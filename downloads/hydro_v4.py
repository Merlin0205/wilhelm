# -*- coding: utf-8 -*-
"""
Hydrostatická kontrola ponorky Wilhelm v4 s náklonem (to, co původní soustředěný model neuměl).

Co počítá
- hmotnost, těžiště, výtlak a střed vztlaku celé sestavy přímo z dílů (web/models/*.stl) a z poloh matic ve v4_stats.json
- dutinu komory a sifonu rozloženou na krychličky 0,4 mm; plyn ji zaplňuje od nejvyššího místa (invazní perkolace se vztlakem),
  takže vyjde, kolik plynu komora udrží při daném podélném sklonu, kde má plyn těžiště a kdy přeteče přes přepad do výstupu
- kolik plynu zůstane po vyfouknutí (nad úrovní okna sifonu, resp. výstupu)
- rovnovážný sklon při daném množství plynu (moment vztlaku plynu, tíhy a odporu při klesání/stoupání), včetně vlivu volné hladiny plynu
- celý cyklus: zbytek plynu -> plnění -> spuštění sifonu; podmínky: prázdná klesá, plná stoupá, sklon v mezích
- Monte Carlo přes výrobní rozptyl (hustota výtisku, hmotnost matic, uvězněný vzduch, kapilarita, poloha středu odporu)
  s napodobením postupu uživatele: přidávat matice M4, dokud neklesne, pak ještě dvě; při velkém sklonu jednu přesunout

Spuštění: ../.venv312/bin/python sim/hydro_v4.py [--mc N] [--plot]
Úhel: kladný = příď nahoru. Souřadnice jako v generátoru (x od zádi k přídi, z nahoru od kýlu), mm, g, cm3.
"""
import numpy as np, trimesh, json, heapq, math, os, sys, argparse
from shapely.geometry import Polygon
import shapely

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); WM = os.path.join(ROOT, "web", "models"); OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
st = json.load(open(os.path.join(ROOT, "v4_stats.json"), encoding="utf-8"))
L = lambda n: trimesh.load(os.path.join(WM, n + ".stl"))
PLA = st["pla_density"]; BAL = st["ballast"]; STEEL = BAL["steel_density"]
VOX = 0.4
HOLE_CAP = 2 * 0.072 * math.cos(math.radians(60)) / 0.001 / 9.81     # kapilární výška dírky Ø2 mm při ustupujícím úhlu 60° [mm vodního sloupce]

# ------------------------------------------------------------------ hmotnostní model
hull = [L("hull_stern"), L("hull_bow")]; void = L("void")
plastic = [(m.volume / 1000, m.center_mass) for m in hull] + [(L(n).volume / 1000, L(n).center_mass) for n in ("tower", "cartridge", "gun", "flak", "lid", "planes", "rudders")]
pn = st["pins"]; vpin = math.pi * (pn["d"] / 2) ** 2 * pn["l"] / 1000
plastic += [(vpin, np.array([pn["x"], s * pn["y"], pn["z"]])) for s in (-1, 1)]
V_PL = sum(v for v, _ in plastic); C_PL = sum(v * np.array(c) for v, c in plastic) / V_PL
V_VOID = void.volume / 1000; C_VOID = np.array(void.center_mass)
# střed odporu při svislém pohybu: těžiště půdorysu trupu (obrys v rovině největší šířky)
_sec = trimesh.util.concatenate(hull).section(plane_origin=(0, 0, 14.6), plane_normal=(0, 0, 1))
_out = max((Polygon(e[:, :2]) for e in _sec.discrete if len(e) > 3), key=lambda p: p.area)
C_DRAG = np.array([_out.centroid.x, 0.0, 16.0])

# pořadí, v jakém uživatel přidává matice M4 na doladění tíhy (z generátoru: střídavě za komoru a na víčko, nakonec záď a příď)
P8 = BAL["m8_pos"]; TRIM_SEQ = [np.array(t["pos"], float) for t in BAL["trim_seq"]]; NOM = BAL["nominal_trim"]
X4 = BAL["m4_slots_x"]                                             # [záď, střed 1, střed 2, příď 1, příď 2]
ADD = 3                                                            # kolik matic přidat poté, co ponorka začne klesat
FAR = {X4[0]: [8, 10], X4[3]: [12, 13]}                             # pozice v pořadí, které leží na zádi a na přídi (pro přesun matice kvůli sklonu)
def mass_model(n_trim=None, extra=(), k_pla=1.0, k_m8=1.0, k_m4=1.0, air=None, m8_fill=None, bub=0.0):
    """vrací M [g], těžiště, výtlak [cm3], střed vztlaku (bez plynu v komoře). extra = další matice M4 [(x_štěrbiny, k)], air = (cm3, poloha)"""
    mass = [(V_PL * PLA * k_pla, C_PL)]; disp = [(V_PL, C_PL), (V_VOID, C_VOID)]
    for x in (m8_fill if m8_fill is not None else BAL["m8_fill_x"]):
        g = BAL["m8_g"] * k_m8; mass.append((g, np.array(P8[str(x)], float))); disp.append((BAL["m8_g"] / STEEL, np.array(P8[str(x)], float)))
    for p in list(TRIM_SEQ[:NOM if n_trim is None else n_trim]) + [TRIM_SEQ[i] for i in extra]:
        g = BAL["m4_g"] * k_m4; mass.append((g, p)); disp.append((BAL["m4_g"] / STEEL, p))
    if air: disp.append(air)
    if bub: disp.append((bub, np.array([166.0, 0.0, 20.0])))          # bublinky CO2 přichycené na trupu při hře
    M = sum(a for a, _ in mass); cg = sum(a * c for a, c in mass) / M
    V = sum(a for a, _ in disp); cb = sum(a * c for a, c in disp) / V
    return M, cg, V, cb

# ------------------------------------------------------------------ dutina komory a sifonu jako krychličky
cu, si, fu = st["cup"], st["siphon"], st["fuel"]
cav = trimesh.boolean.difference([trimesh.boolean.union([L("cup"), L("siphon")]), trimesh.boolean.union([L("cartridge"), L("lid")])])
x0, x1 = cu["x0"] - 1.0, si["xc"] + si["d"] / 2 + 1.0; y0, y1 = -cu["w"] / 2 - 1.0, 13.0; z0, z1 = 0.0, si["exit_z"] + si["d"] / 2 + 1.0
xs = np.arange(x0 + VOX / 2, x1, VOX); ys = np.arange(y0 + VOX / 2, y1, VOX); zs = np.arange(z0 + VOX / 2, z1, VOX)
YY, ZZ = np.meshgrid(ys, zs, indexing="ij")
grid = np.zeros((len(xs), len(ys), len(zs)), bool)
for i, x in enumerate(xs):
    s = cav.section(plane_origin=(x, 0, 0), plane_normal=(1, 0, 0))
    if s is None: continue
    poly = None
    for e in s.discrete:
        if len(e) < 4: continue
        p = Polygon(e[:, 1:3]).buffer(0)
        poly = p if poly is None else poly.symmetric_difference(p)
    if poly is None or poly.is_empty: continue
    grid[i] = shapely.contains_xy(poly, YY, ZZ)
X3, Y3, Z3 = np.meshgrid(xs, ys, zs, indexing="ij")
EXIT = grid & (Y3 > 11.0)                                           # ústí výstupu na boku
V_CAV = grid.sum() * VOX ** 3 / 1000
SEED = (int(np.argmin(abs(xs - fu["x"]))), int(np.argmin(abs(ys - 0.0))), int(np.argmin(abs(zs - 9.0))))
assert grid[SEED], "zdroj plynu (košíček) není v dutině"
WEIR_X = (si["xb"] + si["d"] / 2 - 0.3, si["xc"] - si["d"] / 2 + 0.3)   # stěna mezi větvemi: pod ní je přepad
SILL = (cu["x0"] + cu["l"] + 0.5, si["top_z"] - si["port_h"])         # spodní hrana okna (x, z)
EXITP = (si["xc"], si["exit_z"])

def invade(theta_deg, cap_mm):
    """plyn stoupá od košíčku a plní dutinu shora; vrací objem před přetečením, zbytek po vyfouknutí, těžiště plynu jako funkci objemu, únik víčkem"""
    th = math.radians(theta_deg); c, s_ = math.cos(th), math.sin(th)
    H = Z3 * c + (X3 - cu["cx"]) * s_
    ok = grid.copy()
    ok &= ~((X3 > WEIR_X[0]) & (X3 < WEIR_X[1]) & (Z3 > st["weir_z"] - cap_mm))      # kapilarita: plyn projde až o cap_mm níž
    H = np.where(Z3 < st["plug_t"], H - HOLE_CAP, H)                                  # dírky víčka Ø2 mm: plyn jimi projde až po překonání kapilárního tlaku
    nx, ny, nz = ok.shape; seen = np.zeros(ok.shape, bool); heap = [(-H[SEED], SEED)]; seen[SEED] = True
    order = []; lid_leak = False
    while heap:
        _, (i, j, k) = heapq.heappop(heap)
        if EXIT[i, j, k]: break
        if zs[k] < st["plug_t"]: lid_leak = True; continue             # plyn protlačený dírkami víčka ven (ztráta)
        order.append((i, j, k))
        for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            a, b, d = i + di, j + dj, k + dk
            if 0 <= a < nx and 0 <= b < ny and 0 <= d < nz and ok[a, b, d] and not seen[a, b, d]:
                seen[a, b, d] = True; heapq.heappush(heap, (-H[a, b, d], (a, b, d)))
    o = np.array(order); px, pz = xs[o[:, 0]], zs[o[:, 2]]; hh = pz * c + (px - cu["cx"]) * s_
    dv = VOX ** 3 / 1000; vol = np.arange(1, len(o) + 1) * dv
    cx = np.cumsum(px) / np.arange(1, len(o) + 1); cz = np.cumsum(pz) / np.arange(1, len(o) + 1)
    h_stop = min(SILL[1] * c + (SILL[0] - cu["cx"]) * s_, EXITP[1] * c + (EXITP[0] - cu["cx"]) * s_)
    keep = hh > h_stop; v_ret = keep.sum() * dv
    rx = px[keep].mean() if keep.any() else cu["cx"]; rz = pz[keep].mean() if keep.any() else cu["top"]
    return dict(theta=theta_deg, cap=vol[-1], ret=v_ret, ret_c=(rx, rz), vol=vol, cx=cx, cz=cz, lid_leak=lid_leak)

THETAS = np.arange(-16, 16.1, 2.0)
_cache = {}
def table(cap_mm):
    key = round(cap_mm, 2)
    if key not in _cache: _cache[key] = [invade(t, cap_mm) for t in THETAS]
    return _cache[key]

def gas_state(tab, theta, V):
    """těžiště plynu objemu V při sklonu theta (lineární interpolace mezi vypočtenými sklony); vrací (cx, cz, kapacita, zbytek)"""
    t = float(np.clip(theta, THETAS[0], THETAS[-1])); i = min(int((t - THETAS[0]) // 2), len(THETAS) - 2); w = (t - THETAS[i]) / 2.0
    out = []
    for r in (tab[i], tab[i + 1]):
        k = min(max(int(V / (VOX ** 3 / 1000)) - 1, 0), len(r["vol"]) - 1); out.append((r["cx"][k], r["cz"][k], r["cap"], r["ret"]))
    return tuple((1 - w) * a + w * b for a, b in zip(*out))

def pitch_eq(mm, tab, V, theta0=0.0):
    """rovnovážný sklon [°] pro objem plynu V: vztlak (trup + plyn) a odpor (v těžišti půdorysu) proti tíze v těžišti"""
    M, cg, Vd, cb = mm; th = theta0
    for _ in range(60):
        gx, gz, cap, ret = gas_state(tab, th, V)
        B = Vd + V; D = M - B                                          # odpor při ustáleném klesání (kladný) nebo stoupání (záporný)
        ex = (Vd * cb[0] + V * gx + D * C_DRAG[0]) / M; ez = (Vd * cb[2] + V * gz + D * C_DRAG[2]) / M
        new = math.degrees(math.atan2(ex - cg[0], ez - cg[2])); new = max(-40.0, min(40.0, new))
        if abs(new - th) < 0.02: th = new; break
        th = th + 0.5 * (new - th)                                     # tlumená iterace (volná hladina plynu sklon zesiluje)
    return th

def cycle(mm, cap_mm=3.5):
    """jeden cyklus: od zbytku plynu po spuštění sifonu. Vrací slovník s objemy, tíhou ve vodě a sklonem."""
    tab = table(cap_mm); M, cg, Vd, cb = mm; W0 = M - Vd
    th = 0.0
    for _ in range(20):                                                # stav po vyfouknutí: zbytek plynu a sklon se navzájem ovlivňují
        _, _, cap, ret = gas_state(tab, th, 0.05); th_new = pitch_eq(mm, tab, ret, th)
        if abs(th_new - th) < 0.05: th = th_new; break
        th = th_new
    V_ret, th_ret = ret, th
    V = V_ret; th_max = th_min = th; trig = None
    while V < 12.0:
        V += 0.05; th = pitch_eq(mm, tab, V, th); th_max = max(th_max, th); th_min = min(th_min, th)
        _, _, cap, _ = gas_state(tab, th, V)
        if V >= cap: trig = V; break
    leak = any(r["lid_leak"] for r in tab if th_min - 2 <= r["theta"] <= th_max + 2)
    res = dict(W0=W0, V_ret=V_ret, V_trig=trig, th_ret=th_ret, th_trig=th, th_min=th_min, th_max=th_max, W_ret=W0 - V_ret, W_trig=(W0 - trig) if trig else None, lid_leak=leak)
    res["ok"] = bool(trig and res["W_ret"] > 0.25 and res["W_trig"] < -0.25 and max(abs(th_min), abs(th_max)) <= 15.0 and not leak)
    return res

# ------------------------------------------------------------------ napodobení uživatele
def user_trim(sample, cap_mm):
    """postup z návodu: ubírat/přidávat matice M4, až ponorka (bez paliva, po protřepání) právě klesá, pak přidat tři;
    když klesá nebo stoupá šikmo (střed sklonu přes 6°), jednu matici přesunout na záď nebo na příď"""
    kw = dict(k_pla=sample["k_pla"], k_m8=sample["k_m8"], k_m4=sample["k_m4"], air=sample["air"])
    n = None
    for k in range(0, len(TRIM_SEQ) - ADD + 1):                         # 8 míst uprostřed (štěrbiny za komorou a víčko), pak konce trupu
        mm = mass_model(n_trim=min(k, 8), extra=list(range(8, k)), **kw); c = cycle(mm, cap_mm)
        if mm[0] - mm[2] - c["V_ret"] > 0.0: n = k; break
    if n is None: return None, dict(fail="ani se všemi maticemi M4 neklesne (výtisk je příliš lehký)")
    n = n + ADD if n > 0 else 0                                         # klesá-li už bez matic, nepřidávat nic
    seq_n = min(n, 8); extra = list(range(8, n))                        # nad osm: konce trupu (pozice 8, 9, …)
    mm = mass_model(n_trim=seq_n, extra=extra, **kw); c = cycle(mm, cap_mm)
    for _ in range(2):
        mid = 0.5 * (c["th_min"] + c["th_max"])
        if abs(mid) <= 6.0 or seq_n == 0: break
        cand = [i for i in FAR[X4[0] if mid < 0 else X4[3]] if i not in extra]      # příď dolů (mid < 0) -> matici na záď; příď nahoru -> na příď
        if not cand: break
        seq_n -= 1; extra.append(cand[0])
        mm = mass_model(n_trim=seq_n, extra=extra, **kw); c = cycle(mm, cap_mm)
    def play(): m_ = mass_model(n_trim=seq_n, extra=extra, bub=sample.get("bub", 0.0), **kw); return m_, cycle(m_, cap_mm)     # při hře navíc bublinky na trupu
    mm, c = play(); first_ok = c["ok"]; fixes = 0
    while not c["ok"] and fixes < 3 and c["V_trig"]:                   # rady z návodu: neklesne znovu -> přidat matici; nevyplave -> ubrat matici
        if c["W_ret"] <= 0.25:
            if seq_n < 8: seq_n += 1
            else:
                free = [i for i in range(8, len(TRIM_SEQ)) if i not in extra]
                if not free: break
                extra.append(free[0])
        elif c["W_trig"] >= -0.25:
            if extra: extra.pop()
            elif seq_n > 0: seq_n -= 1
            else: break
        else: break
        fixes += 1; mm, c = play()
    c["n_trim"] = seq_n; c["extra"] = extra; c["first_ok"] = first_ok; c["fixes"] = fixes
    return mm, c

def monte_carlo(N, seed=1):
    rng = np.random.default_rng(seed); res = []; caps = (2.5, 3.5, 4.5)
    for _ in range(N):
        s = dict(k_pla=rng.uniform(0.972, 1.03), bub=rng.uniform(0.0, 0.6),   # hustota výtisku 1,215–1,29 g/cm3 (PLA Basic nebo běžné PLA, podtečení i přetečení při plné výplni); bublinky na trupu až 0,6 cm3
                 k_m8=rng.uniform(0.96, 1.05), k_m4=rng.uniform(0.93, 1.05),
                 air=(rng.uniform(0.0, 0.8), np.array([rng.uniform(60, 280), 0.0, rng.uniform(10, 30)])))
        cap_mm = caps[int(rng.integers(0, 3))]
        mm, c = user_trim(s, cap_mm)
        if mm is None: res.append(dict(ok=False, fail=c["fail"])); continue
        c["fail"] = None if c["ok"] else ("nespustí sifon" if not c["V_trig"] else "neklesá po vyfouknutí" if c["W_ret"] <= 0.25 else "nevyplave" if c["W_trig"] >= -0.25 else "sklon přes 15°" if max(abs(c["th_min"]), abs(c["th_max"])) > 15 else "únik víčkem")
        res.append(c)
    return res

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--mc", type=int, default=0); ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    mm = mass_model(); M, cg, Vd, cb = mm
    print(f"hmotnost {M:.1f} g, výtlak {Vd:.1f} cm3, tíha ve vodě bez plynu {M - Vd:+.2f} g; těžiště x {cg[0]:.2f} z {cg[2]:.2f}; střed vztlaku x {cb[0]:.2f} z {cb[2]:.2f}; BG {cb[2] - cg[2]:.2f} mm; střed odporu x {C_DRAG[0]:.1f}")
    print(f"dutina komory a sifonu {V_CAV:.2f} cm3 ({grid.sum()} krychliček {VOX} mm)")
    tab = table(3.5)
    print(" sklon   udrží plynu   zbude po vyfouknutí   těžiště plynu x")
    for r in tab: print(f" {r['theta']:+5.0f}°   {r['cap']:6.2f} cm3      {r['ret']:5.2f} cm3          {r['cx'][-1]:6.1f}   {'ÚNIK VÍČKEM' if r['lid_leak'] else ''}")
    out = {"mass_g": round(M, 2), "disp_cm3": round(Vd, 2), "bg_mm": round(float(cb[2] - cg[2]), 2), "capacity": [dict(theta=float(r["theta"]), cap=round(float(r["cap"]), 2), ret=round(float(r["ret"]), 2)) for r in tab], "cycles": {}}
    for cap_mm in (2.5, 3.5, 4.5):
        c = cycle(mm, cap_mm); out["cycles"][str(cap_mm)] = {k: (round(float(v), 2) if isinstance(v, (int, float)) and not isinstance(v, bool) else v) for k, v in c.items()}
        print(f"kapilarita {cap_mm} mm: zbytek {c['V_ret']:.2f} cm3 (tíha {c['W_ret']:+.2f} g, sklon {c['th_ret']:+.1f}°) -> spuštění při {c['V_trig']:.2f} cm3 (tíha {c['W_trig']:+.2f} g, sklon {c['th_trig']:+.1f}°); rozsah sklonu {c['th_min']:+.1f}° až {c['th_max']:+.1f}°; {'OK' if c['ok'] else 'NEVYHOVUJE'}")
    if a.mc:
        res = monte_carlo(a.mc); okn = sum(1 for r in res if r["ok"]); ok1 = sum(1 for r in res if r.get("first_ok")); fails = {}
        for r in res:
            if not r["ok"]: fails[r["fail"]] = fails.get(r["fail"], 0) + 1
        good = [r for r in res if r["ok"]]
        print(f"Monte Carlo {a.mc} ponorek: po vyvážení podle návodu funguje napoprvé {ok1} ({100 * ok1 / a.mc:.0f} %), po úpravě podle rad „když něco nejde“ {okn} ({100 * okn / a.mc:.1f} %); selhání {fails}")
        if good:
            print(f"  dolaďovacích matic M4: medián {np.median([r['n_trim'] + len(r['extra']) for r in good]):.0f} (rozsah {min(r['n_trim'] + len(r['extra']) for r in good)}–{max(r['n_trim'] + len(r['extra']) for r in good)}); sklon během cyklu: 90 % ponorek do ±{np.percentile([max(abs(r['th_min']), abs(r['th_max'])) for r in good], 90):.1f}°")
            print(f"  síla při klesání medián {np.median([r['W_ret'] for r in good]):.2f} g (min {min(r['W_ret'] for r in good):.2f}), při stoupání medián {-np.median([r['W_trig'] for r in good]):.2f} g (min {-max(r['W_trig'] for r in good):.2f})")
        out["mc"] = dict(n=a.mc, ok=okn, ok_first=ok1, fails=fails)
    json.dump(out, open(os.path.join(OUT, "hydro_v4.json"), "w"), indent=1, ensure_ascii=False)
