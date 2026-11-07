# -*- coding: utf-8 -*-
"""
Wilhelm v2 – hračka ponorky U-47 (typ VIIB) na prášek do pečiva, měřítko 1:200 (332,5 mm),
dělený trup pro Bambu Lab P1S (podložka 256 × 256 mm).

Princip: komora pod kýlem se plní CO2, při plné komoře ji plynový sifon vypustí celou.
Trup: plný materiál (100 % výplň) + uzavřená dutina v palubní nástavbě jako plovák,
dva díly spojené kolíky. Věž, dělo, zátka, kormidla a kolíky zvlášť.

Spuštění:  python u47_v2.py   (numpy, scipy, trimesh, manifold3d, shapely, networkx, matplotlib)
Výstup:    stl_v2/*.stl, web/models/*.stl (pro prohlížeč), v2_stats.json, v2_sections.json
Souřadnice: x od zádi (0) k přídi (L), y doprava (pravobok +), z nahoru, kýl z=0. mm.
"""
import numpy as np, trimesh, os, json
from scipy.interpolate import PchipInterpolator
from trimesh.creation import box, cylinder, extrude_polygon
from shapely.geometry import Polygon as SPoly, Point
from shapely.ops import unary_union

# ---------------------------------------------------------------- parametry
SCALE = 1000.0 / 200.0        # mm na metr (1:200) = 5
L_REAL = 66.5                 # m (u47.org)
L = L_REAL * SCALE            # 332,5 mm
PLA = 1.24                    # g/cm3 (PETG 1,27: o ~3 g těžší, doladí se šroubky)
TARGET_NET_G = -1.5           # hmotnost minus vztlak bez šroubků a prášku (záporné = plave)
SPLIT_X = 215.0               # dělicí rovina trupu (záďová část 215 mm, příďová 117,5 mm)

CUP_W, CUP_L, CUP_TOP = 12.0, 60.0, 12.0    # komora na prášek
PLUG_T = 2.4                                # zátka s mřížkou
SIPHON_D = 6.0                              # Ø kanálku: Eötvös > 3,4, bubliny stoupají samy, kapilární tlak ~1,7 mm H2O
BEND_Z = 5.5                                # spodek ohybu sifonu (> PLUG_T + kapilární rezerva)
WALL_SIDE, WALL_DECK = 2.2, 2.0             # stěny dutiny plováku
PIN_D, PIN_L = 4.0, 24.0                    # spojovací kolíky

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "stl_v2"); os.makedirs(OUT, exist_ok=True)
WEB = os.path.join(HERE, "web", "models"); os.makedirs(WEB, exist_ok=True)

# ------------------------------------------------- linie trupu (m skutečnosti, x od zádi)
xs_k  = [0, 2, 4, 6, 8, 11, 13, 50, 54, 58, 61, 64, 66.0, 66.5]
zk    = [4.1, 3.2, 2.4, 1.7, 1.0, 0.3, 0, 0, 0.3, 1.2, 2.3, 3.4, 4.7, 5.6]
xs_w  = [0, 2, 5, 8, 14, 20, 26, 40, 46, 52, 58, 62, 65, 66.5]
wm    = [0.05, 0.3, 0.75, 1.2, 1.65, 2.25, 2.95, 3.1, 2.65, 1.95, 1.25, 0.7, 0.28, 0.08]
xs_c  = [0, 2, 4, 10, 18, 26, 44, 52, 58, 63, 66.5]
wc    = [0.15, 0.45, 0.75, 1.15, 1.55, 1.75, 1.75, 1.5, 1.1, 0.6, 0.12]
ZD = 6.5                                     # paluba rovná (tisk bez podpěr)
f_zk, f_wm, f_wc = PchipInterpolator(xs_k, zk), PchipInterpolator(xs_w, wm), PchipInterpolator(xs_c, wc)
def f_zd(x): return ZD

def section(x_m, keel_flat_half=0.0):
    zk_, zd_ = float(f_zk(x_m)), ZD
    wm_, wc_ = float(f_wm(x_m)), float(f_wc(x_m))
    wm_ = max(wm_, wc_ + 0.05)
    h = zd_ - zk_; zm = zk_ + 0.45 * h; zs = zk_ + 0.74 * h
    n = 1.75 if 15 < x_m < 52 else 1.5
    K = 32
    t = np.linspace(0, 1, K); lv = zd_ - (zd_ - zk_) * (t ** 0.8)
    pts = []
    for z in lv:
        if z >= zs: w = wc_
        elif z >= zm:
            tt = (z - zm) / (zs - zm); w = wc_ + (wm_ - wc_) * np.sqrt(max(0.0, 1 - tt * tt))
        else:
            tt = (zm - z) / (zm - zk_); w = wm_ * (max(0.0, 1 - tt ** n)) ** (1.0 / n)
        if keel_flat_half > 0:
            zb = 0.75
            wb = keel_flat_half if z - zk_ < zb else keel_flat_half * max(0.0, 1 - (z - zk_ - zb) / (1.4 * zb))
            w = max(w, wb)
        pts.append((max(w, 0.015), z))
    pts = np.array(pts); left = pts[::-1].copy(); left[:, 0] *= -1
    return np.vstack([pts, left]) * SCALE

def loft(stations_m, flat_fn):
    secs = [section(x, flat_fn(x)) for x in stations_m]; M = len(secs[0])
    V = np.vstack([np.column_stack([np.full(M, x * SCALE), s[:, 0], s[:, 1]]) for x, s in zip(stations_m, secs)])
    F = []
    for i in range(len(secs) - 1):
        a, b = i * M, (i + 1) * M
        for j in range(M):
            j2 = (j + 1) % M; F += [[a + j, b + j, b + j2], [a + j, b + j2, a + j2]]
    c0 = len(V); V = np.vstack([V, V[:M].mean(0)]); F += [[c0, (j + 1) % M, j] for j in range(M)]
    base = (len(secs) - 1) * M; c1 = len(V); V = np.vstack([V, V[base:base + M].mean(0)])
    F += [[c1, base + j, base + (j + 1) % M] for j in range(M)]
    m = trimesh.Trimesh(V, np.array(F), process=True); m.fix_normals(); return m

stations = np.concatenate([np.linspace(0, 3, 12), np.linspace(3.3, 63, 220), np.linspace(63.2, L_REAL, 16)])
hull0 = loft(stations, lambda x: 0.0)
cup_cx = float(hull0.center_mass[0])
tower_cx = (L_REAL - 33.8) * SCALE
print(f"objemové těžiště x = {cup_cx:.1f} mm, věž x = {tower_cx:.1f} mm, L = {L:.1f} mm")
flat_half_m = (CUP_W / 2 + 2.0) / SCALE
def flat_fn(x_m):
    d = abs(x_m * SCALE - cup_cx) - CUP_L / 2
    return flat_half_m if d <= 0 else flat_half_m * max(0.0, 1 - d / 14.0)
hull = loft(stations, flat_fn)
under = trimesh.intersections.slice_mesh_plane(hull, plane_origin=(0, 0, 4.7 * SCALE), plane_normal=(0, 0, -1), cap=True)
print(f"trup vodotěsný {hull.is_watertight}, objem {hull.volume/1000:.1f} cm3, pod vodoryskou {under.volume/1000:.1f} cm3 (753 t -> {753e6/200**3:.1f})")

def B(sx, sy, sz, cx, cy, cz): b = box((sx, sy, sz)); b.apply_translation((cx, cy, cz)); return b
def Cz(d, h, cx, cy, cz): c = cylinder(radius=d / 2, height=h, sections=40); c.apply_translation((cx, cy, cz)); return c
def Cx(d, h, cx, cy, cz):
    c = cylinder(radius=d / 2, height=h, sections=40); c.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, (0, 1, 0))); c.apply_translation((cx, cy, cz)); return c
def Cy(d, h, cx, cy, cz):
    c = cylinder(radius=d / 2, height=h, sections=40); c.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, (1, 0, 0))); c.apply_translation((cx, cy, cz)); return c
def stadium(lx, wy):
    return unary_union([Point(-(lx - wy) / 2, 0).buffer(wy / 2, 32), Point((lx - wy) / 2, 0).buffer(wy / 2, 32),
                        SPoly([(-(lx - wy) / 2, -wy / 2), ((lx - wy) / 2, -wy / 2), ((lx - wy) / 2, wy / 2), (-(lx - wy) / 2, wy / 2)])])
U = trimesh.boolean.union; D = lambda a, b: trimesh.boolean.difference([a, b])

# -------- kormidla (dvě svislá, typ VIIB) – součást záďového dílu
rud = [B(2.0 * SCALE, 1.6, 2.6 * SCALE, 2.6 * SCALE, s * 1.1 * SCALE, float(f_zk(2.6)) * SCALE + 0.3) for s in (-1, 1)]
solid = U([hull] + rud)

# -------- otevřené dutiny
cup = B(CUP_L, CUP_W, CUP_TOP + 1, cup_cx, 0, (CUP_TOP + 1) / 2 - 1)
notch = B(4.0, CUP_W + 4.0, 2.4, cup_cx + CUP_L / 2 + 1.5, 0, 0.0)                   # zářez pro nehet (přední konec)
xa = cup_cx - CUP_L / 2 + 4.0; xb = xa - 9.5; xc = xb - 9.5
z_top = CUP_TOP + 6.0; z_exit = 15.0; r = SIPHON_D / 2   # horní větev 3 mm nad stropem komory; r = SIPHON_D / 2
siphon = U([
    Cz(SIPHON_D, z_top - CUP_TOP + 1 + r, xa, 0, (CUP_TOP - 1 + z_top + r) / 2),
    Cx(SIPHON_D, xa - xb + SIPHON_D, (xa + xb) / 2, 0, z_top),
    Cz(SIPHON_D, z_top - BEND_Z + SIPHON_D, xb, 0, (z_top + BEND_Z) / 2),
    Cx(SIPHON_D, xb - xc + SIPHON_D, (xb + xc) / 2, 0, BEND_Z),
    Cz(SIPHON_D, z_exit - BEND_Z + SIPHON_D, xc, 0, (z_exit + BEND_Z) / 2),
    Cy(SIPHON_D, 28.0, xc, 14.0, z_exit),
])
# šroubky M4 (doladění zátěže a trimu): Ø3,3 × 12 mm, dva vzadu, dva vpředu
screw_x = [xc - 14.0, xc - 42.0, cup_cx + CUP_L / 2 + 14.0, SPLIT_X + 30.0]
screws = [Cz(3.3, 13.0, x, 0, 5.5) for x in screw_x]
deck_z = ZD * SCALE
tower_recess = extrude_polygon(stadium(6.4 * SCALE + 0.4, 1.9 * SCALE + 0.4), 0.8); tower_recess.apply_translation((tower_cx, 0, deck_z - 0.4))
gun_x = (L_REAL - 26.0) * SCALE; flak_x = (L_REAL - 40.0) * SCALE
gun_hole = Cz(2.3, 5.0, gun_x, 0, deck_z - 2.0); flak_hole = Cz(1.8, 4.0, flak_x, 0, deck_z - 1.5)
slot_bow = B(7.2, 40, 1.9, 60.5 * SCALE, 0, 3.3 * SCALE)
slot_stern = B(8.2, 40, 1.9, 4.0 * SCALE, 0, 3.4 * SCALE)
pin_holes = [Cx(PIN_D + 0.15, PIN_L + 2.0, SPLIT_X, s * 4.5, 12.0) for s in (-1, 1)]
hatches = [Cz(4.0, 1.0, x, 0, deck_z) for x in (12.0 * SCALE, 17.5 * SCALE, 55.0 * SCALE)]
# zaplavovací otvory nástavby: řada obdélníkových štěrbin po obou bocích
slots = []
for x in np.arange(7.0, 63.5, 1.4):
    if abs(x * SCALE - tower_cx) < 18: pass
    y = float(f_wc(x)) * SCALE
    for s in (-1, 1):
        slots.append(B(4.0, 1.8, 1.3, x * SCALE, s * y, deck_z - 2.6))
open_cav = U([cup, notch, siphon, tower_recess, gun_hole, flak_hole, slot_bow, slot_stern] + screws + pin_holes + hatches + slots)
A = D(solid, open_cav)
V_disp = A.volume / 1000.0

# -------- uzavřená dutina plováku v nástavbě
def make_void(z_lo):
    segs = []
    for x in stations:
        wc_ = float(f_wc(x)) * SCALE - WALL_SIDE; zd_ = deck_z - WALL_DECK
        if wc_ > 1.5 and zd_ - z_lo > 1.0: segs.append((x * SCALE, wc_, zd_))
    boxes = []
    for i in range(len(segs) - 1):
        x0, x1 = segs[i][0], segs[i + 1][0]; w = min(segs[i][1], segs[i + 1][1]); zt = min(segs[i][2], segs[i + 1][2])
        boxes.append(B(x1 - x0 + 0.01, 2 * w, zt - z_lo, (x0 + x1) / 2, 0, (zt + z_lo) / 2))
    v = U(boxes)
    keep = U([
        B(xa - xc + 2 * SIPHON_D + 6.0, 60, 24.0, (xa + xc) / 2, 0, 12.0),   # sifon
        B(6.0, 60, 60, SPLIT_X, 0, 20),                                     # stěny u dělicí roviny (3 mm na každé straně)
        Cz(8.0, 40, gun_x, 0, 20), Cz(7.0, 40, flak_x, 0, 20),              # pod dělem a flakem
    ])
    return D(v, keep)
def net_for(z_lo):
    v = make_void(z_lo); H = D(A, v); return H, v, H.volume / 1000.0 * PLA - V_disp

# -------- věž
tower_h = (9.5 - 6.5) * SCALE; tw_l, tw_w = 6.4 * SCALE, 1.9 * SCALE
tower = extrude_polygon(stadium(tw_l, tw_w), tower_h)
well = extrude_polygon(stadium(tw_l - 3.2, tw_w - 3.2), 3.0); well.apply_translation((0, 0, tower_h - 3.0 + 0.01))
tower = D(tower, well)
scopes = [Cz(2.2, 9.0, -4.0, 1.2, tower_h - 3.0 + 4.5), Cz(1.8, 7.0, -0.5, -1.4, tower_h - 3.0 + 3.5), Cz(1.6, 3.5, 9.0, 0, tower_h - 3.0 + 1.75)]
tower = U([tower] + scopes)
# dělo 8,8 cm (tiskne se naležato), čep Ø2 do paluby
gun = U([Cz(4.0, 2.5, 0, 0, 1.25), Cz(2.0, 3.0, 0, 0, -1.5 + 0.01), B(5.0, 3.0, 3.0, 0, 0, 4.0), Cx(1.6, 16.0, 9.0, 0, 4.6), B(2.6, 2.0, 1.6, -3.0, 0, 4.6)])
flak = U([Cz(3.0, 1.5, 0, 0, 0.75), Cz(1.5, 2.5, 0, 0, -1.25 + 0.01), Cz(1.4, 4.0, 0, 0, 3.5), Cx(1.0, 7.0, 3.0, 0, 5.4)])
pins = U([Cx(PIN_D - 0.1, PIN_L, 0, s * 6, PIN_D / 2) for s in (-1, 1)])
tower_net = tower.volume / 1000 * (PLA - 1.0); gun_net = (gun.volume + flak.volume) / 1000 * (PLA - 1.0)
plug_net = CUP_L * CUP_W * PLUG_T / 1000 * (PLA - 1.0) * 0.85
pins_net = pins.volume / 1000 * PLA        # kolíky jsou uvnitř trupu, vodu navíc nevytlačí (otvory se počítají jako dutina)

target = TARGET_NET_G - tower_net - gun_net - plug_net - pins_net
z_a, z_b = 22.0, 27.0
H_a, v_a, n_a = net_for(z_a); H_b, v_b, n_b = net_for(z_b)
for _ in range(7):
    z_new = float(np.clip(z_a + (target - n_a) * (z_b - z_a) / (n_b - n_a), 14.0, deck_z - WALL_DECK - 2.0))
    H_n, v_n, n_n = net_for(z_new)
    z_a, n_a, H_a, v_a = z_b, n_b, H_b, v_b; z_b, n_b, H_b, v_b = z_new, n_n, H_n, v_n
    if abs(n_n - target) < 0.05: break
hull_final, void, z_void = H_b, v_b, z_b
net_total = n_b + tower_net + gun_net + plug_net + pins_net
cup_cap = CUP_L * CUP_W * (CUP_TOP - BEND_Z) / 1000
print(f"dutina od z={z_void:.1f} mm, objem {void.volume/1000:.1f} cm3; vytlačeno {V_disp:.1f} cm3; PLA trup {hull_final.volume/1000:.1f} cm3 = {hull_final.volume/1000*PLA:.0f} g")
print(f"bilance bez šroubků a prášku {net_total:+.2f} g; kapacita plynu {cup_cap:.1f} cm3")
cm, cb = hull_final.center_mass, A.center_mass
print(f"těžiště z={cm[2]:.1f}, střed vztlaku z={cb[2]:.1f} mm")

# -------- rozdělení trupu
big = 1000
stern = trimesh.boolean.intersection([hull_final, B(SPLIT_X + 1, big, big, (SPLIT_X + 1) / 2 - 1, 0, 0)])
bow = trimesh.boolean.intersection([hull_final, B(big, big, big, SPLIT_X + big / 2, 0, 0)])
print("díly vodotěsné:", stern.is_watertight, bow.is_watertight, "délky:", round(stern.extents[0], 1), round(bow.extents[0], 1))

def export(mesh, name, flip=False):
    m = mesh.copy()
    if flip:
        m.apply_transform(trimesh.transformations.rotation_matrix(np.pi, (1, 0, 0)))
    m.apply_translation(-m.bounds[0]); m.export(os.path.join(OUT, name + ".stl")); return m
export(stern, "hull_stern_print", flip=True); export(bow, "hull_bow_print", flip=True)
export(tower, "tower"); export(gun, "gun_88"); export(flak, "flak_20"); export(pins, "pins")
lids = {}
for name, wfit in (("lid_tight", CUP_W + 0.1), ("lid_normal", CUP_W - 0.05), ("lid_loose", CUP_W - 0.25)):
    lid = B(CUP_L - 0.4, wfit, PLUG_T, 0, 0, PLUG_T / 2)
    holes = [Cz(2.0, PLUG_T + 2, x, y, PLUG_T / 2) for x in np.linspace(-CUP_L / 2 + 5, CUP_L / 2 - 5, 10) for y in (-3.0, 3.0)]
    ribs = [B(8, 0.6, PLUG_T - 0.8, x, s * (wfit / 2), PLUG_T / 2) for x in (-20, 0, 20) for s in (-1, 1)]
    lid = D(U([lid] + ribs), U(holes)); lid.export(os.path.join(OUT, name + ".stl")); lids[name] = lid
bow_plane = B(7.0, 30.0, 1.6, 0, 0, 0.8); stern_plane = B(8.0, 28.0, 1.6, 12, 0, 0.8)
U([bow_plane, stern_plane]).export(os.path.join(OUT, "planes.stl"))

# -------- zkušební kostka komory se sifonem (tisk ~30 min): stejná komora, stejný sifon, stejná zátka
tb_x0 = xc - 7.0; tb_x1 = cup_cx + CUP_L / 2 + 6.0
test_block = B(tb_x1 - tb_x0, CUP_W + 8.0, CUP_TOP + 12.0, (tb_x0 + tb_x1) / 2, 0, (CUP_TOP + 12.0) / 2)
test_block = D(test_block, U([cup, notch, siphon]))
tb = test_block.copy(); tb.apply_transform(trimesh.transformations.rotation_matrix(np.pi, (1, 0, 0))); tb.apply_translation(-tb.bounds[0]); tb.export(os.path.join(OUT, "test_block_print.stl"))
print("zkušební kostka:", [round(v, 1) for v in test_block.extents], "mm, PLA", round(test_block.volume / 1000 * PLA, 1), "g")

# -------- modely pro web (v poloze sestavení, souřadnice modelu)
for mesh, name in ((stern, "hull_stern"), (bow, "hull_bow"), (void, "void"), (cup, "cup"), (siphon, "siphon")):
    mesh.export(os.path.join(WEB, name + ".stl"))
tower.copy().apply_translation((tower_cx, 0, deck_z)).export(os.path.join(WEB, "tower.stl"))
gun.copy().apply_translation((gun_x, 0, deck_z)).export(os.path.join(WEB, "gun.stl"))
flak.copy().apply_translation((flak_x, 0, deck_z)).export(os.path.join(WEB, "flak.stl"))
lids["lid_normal"].copy().apply_translation((cup_cx, 0, 0)).export(os.path.join(WEB, "lid.stl"))
bp = B(7.0, 30.0, 1.6, 60.5 * SCALE, 0, 3.3 * SCALE); sp = B(8.0, 28.0, 1.6, 4.0 * SCALE, 0, 3.4 * SCALE)
U([bp, sp]).export(os.path.join(WEB, "planes.stl"))

# -------- řezy pro schéma
def sec_paths(m, origin, normal, ax):
    s = m.section(plane_origin=origin, plane_normal=normal); out = []
    for e in s.discrete:
        pts = e[:, ax]; keep = [pts[0]]
        for p in pts[1:]:
            if np.linalg.norm(p - keep[-1]) > 0.2: keep.append(p)
        out.append([[round(float(p[0]), 2), round(float(p[1]), 2)] for p in keep])
    return out
json.dump({"long": sec_paths(hull_final, (0, 0, 0), (0, 1, 0), [0, 2]),
           "cross_cup": sec_paths(hull_final, (cup_cx, 0, 0), (1, 0, 0), [1, 2]),
           "cross_siphon": sec_paths(hull_final, (xb, 0, 0), (1, 0, 0), [1, 2])}, open(os.path.join(HERE, "v2_sections.json"), "w"))

g = lambda m: round(m.volume / 1000 * PLA, 1)
stats = {
    "scale": 200, "length_mm": round(L, 1), "beam_mm": round(2 * 3.1 * SCALE, 1), "height_mm": round(9.5 * SCALE, 1), "deck_mm": round(deck_z, 1),
    "waterline_mm": round(4.7 * SCALE, 1), "displaced_cm3": round(V_disp, 1), "void_cm3": round(void.volume / 1000, 1), "void_z_lo_mm": round(z_void, 1),
    "mass_g": {"hull_stern": g(stern), "hull_bow": g(bow), "tower": g(tower), "gun": g(gun), "flak": g(flak), "lid": round(CUP_L * CUP_W * PLUG_T / 1000 * PLA * 0.85, 1),
               "planes": g(U([bow_plane, stern_plane])), "pins": g(pins)},
    "net_g": round(net_total, 2), "gas_capacity_cm3": round(cup_cap, 1), "cg_z_mm": round(cm[2], 1), "cb_z_mm": round(cb[2], 1),
    "cup": {"w": CUP_W, "l": CUP_L, "top": CUP_TOP, "cx": round(cup_cx, 1)}, "siphon": {"d": SIPHON_D, "bend_z": BEND_Z, "exit_z": z_exit, "xa": round(xa, 1), "xb": round(xb, 1), "xc": round(xc, 1)},
    "split_x": SPLIT_X, "tower_cx": round(tower_cx, 1), "gun_x": round(gun_x, 1), "flak_x": round(flak_x, 1), "screw_x": [round(x, 1) for x in screw_x],
    "parts_extents_mm": {"hull_stern": [round(v, 1) for v in stern.extents], "hull_bow": [round(v, 1) for v in bow.extents], "tower": [round(v, 1) for v in tower.extents]},
}
stats["mass_g"]["total"] = round(sum(stats["mass_g"].values()), 1)
json.dump(stats, open(os.path.join(HERE, "v2_stats.json"), "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats["mass_g"], ensure_ascii=False), "| díly:", stats["parts_extents_mm"])
print("hotovo:", OUT)
