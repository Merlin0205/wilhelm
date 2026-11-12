# -*- coding: utf-8 -*-
"""
Wilhelm v4 – ponorka na prášek do pečiva: U-47 (typ VIIB) 1:200, 332,5 mm, dělený trup pro Bambu Lab P1S.

Co je nového proti v3 (revize 12. 11. 2026, důvody v CLAUDE.md):
- ocelová zátěž v kýlu (matice M8 a M4 na hraně ve štěrbinách skříňového kýlu) -> těžiště 3–4 mm pod středem vztlaku
- těžiště a střed vztlaku nad sebou (řeší se délkou vzduchové dutiny), dutina kopíruje tvar trupu (stěna 2 mm)
- komora krátká a vysoká (30 × 12 × 19 mm) u středu vztlaku: malá volná hladina plynu, necitlivá na náklon
- sifon: vstup oknem u stropu komory, sestupná větev, přepad 8 mm nad kýlem, vzestupná větev, výstup na pravobok
- západky kazety přitlačí věž k palubě (těsní šachtu), košíček je z boku otevřený (jde do něj vložit kousek tablety)
- hloubková kormidla ze dvou dílů (dřív nešla prostrčit trupem)
(původní text v3 následuje)

Co je nového proti v2 (podle simulace a rešerše):
- sifon Ø 9 mm (Eötvös 11: bubliny stoupají samy, ohyb se po vyfouknutí vždy znovu zalije), komora 14 mm, ohyb 7 mm
- komora posunutá 25 mm před střed vztlaku: plyn zvedá příď při stoupání, šroubky vpředu sklání příď při klesání
  -> s natočenými hloubkovými kormidly klouže vpřed na obou fázích (podvodní kluzák), kormidla zatáčí
- hloubková kormidla otočná na hřídeli Ø 4 se západkami 0/±10/±20°, kormidla otočná na čepech Ø 3
- reliéfní nápis WILHELM na bocích trupu, U 47 a funící býk na věži
- dělicí rovina posunuta na 240 mm (komora celá v záďovém dílu)

Spuštění: ../.venv312/bin/python u47_v4.py  -> stl_v4/, web/models/, v4_stats.json, v4_sections.json
Souřadnice: x od zádi (0) k přídi (L), y pravobok +, z nahoru, kýl z=0, mm.
"""
import numpy as np, trimesh, os, json, math
from scipy.interpolate import PchipInterpolator
from trimesh.creation import box, cylinder, extrude_polygon
from shapely.geometry import Polygon as SPoly, Point
from shapely.ops import unary_union
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

SCALE = 1000.0 / 200.0; L_REAL = 66.5; L = L_REAL * SCALE
PLA = 1.25                              # hustota výtisku: Bambu PLA Basic 1,26 g/cm3 při ~99 % zaplnění (rozptyl se doladí maticemi M4)
TARGET_W_G = 2.0                        # tíha ve vodě s jmenovitou zátěží a komorou bez plynu (kladná = klesá); půlka zdvihu plynu
SPLIT_X = 262.0
CH_L, CH_W = 30.0, 12.0; CH_TOP_R, CH_TOP_F = 18.0, 19.5   # komora: strop stoupá k oknu sifonu
CUP_W, CUP_L, CUP_TOP = CH_W, CH_L, CH_TOP_F
PLUG_T = 2.4
KEEL_H = 7.4                            # výška skříňového kýlu se svislými boky (v něm leží matice naplocho)
SIPHON_D = 9.0; WEIR_Z = 9.0; KEEL_T = 2.0; EXIT_Z = 20.0; PORT_H = 2.0
BEND_Z = WEIR_Z
WALL = 1.8
M8 = dict(af=13.0, t=6.5, e=15.0, g=5.2); M4 = dict(af=7.0, t=3.2, e=7.7, g=0.81); STEEL = 7.85
M8_X = [74.5, 91.5, 108.5, 125.5, 226.0]                    # matice M8: hlavní zátěž (všechny osazené)
M4_X = [62.0, 145.5, 154.2, 240.0, 250.0]                   # matice M4, každá štěrbina na dvě nad sebou: uprostřed doladění tíhy, na koncích doladění sklonu
M8_FILL = list(M8_X); NOMINAL_TRIM = 4                      # jmenovité osazení: 5× M8 a první 4 matice M4 z pořadí TRIM_SEQ (dvě ve štěrbině za komorou, dvě na víčku)
LID_TUBE_X = 8.25                                           # trubička na matice M4 na víčku komory (od středu komory dopředu)      # jmenovité osazení (ostatní štěrbiny jsou rezerva na doladění)
PIN_D, PIN_L = 4.0, 24.0
PLANE_SHAFT_D = 4.0; PLANE_SPAN = 34.0; PLANE_CHORD = 9.0; PLANE_T = 1.6
RUDDER_PIN_D = 3.0
FONT = "/System/Library/Fonts/Supplemental/DIN Condensed Bold.ttf"

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "stl_v4"); os.makedirs(OUT, exist_ok=True)
WEB = os.path.join(HERE, "web", "models"); os.makedirs(WEB, exist_ok=True)

# ---------------- linie trupu (m skutečnosti, x od zádi)
xs_k = [0, 2, 4, 6, 8, 11, 13, 50, 54, 58, 61, 64, 66.0, 66.5]; zk = [4.1, 3.2, 2.4, 1.7, 1.0, 0.3, 0, 0, 0.3, 1.2, 2.3, 3.4, 4.7, 5.6]
xs_w = [0, 2, 5, 8, 14, 20, 26, 40, 46, 52, 58, 62, 65, 66.5]; wm = [0.05, 0.3, 0.75, 1.2, 1.65, 2.25, 2.95, 3.1, 2.65, 1.95, 1.25, 0.7, 0.28, 0.08]
xs_c = [0, 2, 4, 10, 18, 26, 44, 52, 58, 63, 66.5]; wc = [0.15, 0.45, 0.75, 1.15, 1.55, 1.75, 1.75, 1.5, 1.1, 0.6, 0.12]
ZD = 6.5
f_zk, f_wm, f_wc = PchipInterpolator(xs_k, zk), PchipInterpolator(xs_w, wm), PchipInterpolator(xs_c, wc)

def half_width(x_m, z_mm, keel_flat_half=0.0):
    """půlšířka trupu (mm) ve stanici x (m) ve výšce z (mm nad kýlem) – stejný předpis jako řez"""
    zk_, zd_ = float(f_zk(x_m)), ZD; wm_, wc_ = float(f_wm(x_m)), float(f_wc(x_m)); wm_ = max(wm_, wc_ + 0.05)
    h = zd_ - zk_; zm = zk_ + 0.45 * h; zs = zk_ + 0.74 * h; n = 1.75 if 15 < x_m < 52 else 1.5
    z = z_mm / SCALE
    if z >= zs: w = wc_
    elif z >= zm: tt = (z - zm) / (zs - zm); w = wc_ + (wm_ - wc_) * math.sqrt(max(0.0, 1 - tt * tt))
    else: tt = (zm - z) / max(zm - zk_, 1e-6); w = wm_ * (max(0.0, 1 - tt ** n)) ** (1.0 / n)
    if keel_flat_half > 0:
        zb = KEEL_H / SCALE; wb = keel_flat_half if z - zk_ < zb else keel_flat_half * max(0.0, 1 - (z - zk_ - zb) / 0.8); w = max(w, wb)
    return max(w, 0.015) * SCALE

def section(x_m, keel_flat_half=0.0):
    zk_, zd_ = float(f_zk(x_m)), ZD; K = 32; t = np.linspace(0, 1, K); lv = zd_ - (zd_ - zk_) * (t ** 0.8)
    pts = np.array([(half_width(x_m, z * SCALE, keel_flat_half), z * SCALE) for z in lv])
    left = pts[::-1].copy(); left[:, 0] *= -1
    return np.vstack([pts, left])

def loft(stations_m, flat_fn):
    secs = [section(x, flat_fn(x)) for x in stations_m]; M = len(secs[0])
    V = np.vstack([np.column_stack([np.full(M, x * SCALE), s[:, 0], s[:, 1]]) for x, s in zip(stations_m, secs)])
    F = []
    for i in range(len(secs) - 1):
        a, b = i * M, (i + 1) * M
        for j in range(M):
            j2 = (j + 1) % M; F += [[a + j, b + j, b + j2], [a + j, b + j2, a + j2]]
    c0 = len(V); V = np.vstack([V, V[:M].mean(0)]); F += [[c0, (j + 1) % M, j] for j in range(M)]
    base = (len(secs) - 1) * M; c1 = len(V); V = np.vstack([V, V[base:base + M].mean(0)]); F += [[c1, base + j, base + (j + 1) % M] for j in range(M)]
    m = trimesh.Trimesh(V, np.array(F), process=True); m.fix_normals(); return m

stations = np.concatenate([np.linspace(0, 3, 12), np.linspace(3.3, 63, 220), np.linspace(63.2, L_REAL, 16)])
hull0 = loft(stations, lambda x: 0.0)
cb_x = float(hull0.center_mass[0])
tower_cx = (L_REAL - 33.8) * SCALE

FUEL_L, FUEL_W = 14.0, 7.2                 # průřez čepele kazety
FUEL_X = tower_cx + 5.0                     # šachta celá pod věží
ch_x0 = FUEL_X - (FUEL_L + 0.4) / 2 - 1.2; ch_x1 = ch_x0 + CH_L; cup_cx = (ch_x0 + ch_x1) / 2
r = SIPHON_D / 2
xa = ch_x1; xb = ch_x1 + 2.0 + r; xc = xb + SIPHON_D + 3.0
def _ramp(x, a_, b_, w, taper=14.0): d = max(a_ - x, x - b_, 0.0); return w * max(0.0, 1 - d / taper)
def flat_fn(x_m):
    x = x_m * SCALE
    return max(_ramp(x, 57.0, 250.0, 5.5), _ramp(x, 66.0, 250.0, 8.0)) / SCALE                              # skříňový kýl 16 mm (na zádi 11 mm): komora a matice naplocho
hull = loft(stations, flat_fn)
under = trimesh.intersections.slice_mesh_plane(hull, plane_origin=(0, 0, 4.7 * SCALE), plane_normal=(0, 0, -1), cap=True)
print(f"trup vodotěsný {hull.is_watertight}, objem {hull.volume/1000:.1f} cm3, pod vodoryskou {under.volume/1000:.1f} cm3")

def B(sx, sy, sz, cx, cy, cz): b = box((sx, sy, sz)); b.apply_translation((cx, cy, cz)); return b
def Cz(d, h, cx, cy, cz): c = cylinder(radius=d / 2, height=h, sections=48); c.apply_translation((cx, cy, cz)); return c
def Cx(d, h, cx, cy, cz):
    c = cylinder(radius=d / 2, height=h, sections=48); c.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, (0, 1, 0))); c.apply_translation((cx, cy, cz)); return c
def Cy(d, h, cx, cy, cz):
    c = cylinder(radius=d / 2, height=h, sections=48); c.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, (1, 0, 0))); c.apply_translation((cx, cy, cz)); return c
def Sph(d, cx, cy, cz): s = trimesh.creation.icosphere(subdivisions=2, radius=d / 2); s.apply_translation((cx, cy, cz)); return s
def stadium(lx, wy):
    return unary_union([Point(-(lx - wy) / 2, 0).buffer(wy / 2, 32), Point((lx - wy) / 2, 0).buffer(wy / 2, 32),
                        SPoly([(-(lx - wy) / 2, -wy / 2), ((lx - wy) / 2, -wy / 2), ((lx - wy) / 2, wy / 2), (-(lx - wy) / 2, wy / 2)])])
U = trimesh.boolean.union; D = lambda a, b: trimesh.boolean.difference([a, b])

# ---------------- text a znak jako polygony (x, z) -> vytažení podél y
def text_polys(txt, height_mm, font=FONT):
    fp = FontProperties(fname=font) if os.path.exists(font) else FontProperties(family="DejaVu Sans", weight="bold")
    tp = TextPath((0, 0), txt, size=height_mm, prop=fp); polys = [SPoly(p) for p in tp.to_polygons() if len(p) >= 3]
    geom = unary_union([p for p in polys if p.is_valid]) if polys else None
    # díry (vnitřky písmen) – TextPath je vrací jako samostatné polygony; unary_union je spojí, proto použijeme symetrický rozdíl
    from shapely.geometry import MultiPolygon
    res = None
    for p in polys:
        if not p.is_valid: p = p.buffer(0)
        res = p if res is None else res.symmetric_difference(p)
    return res
def bull_poly(width_mm):
    from svgpathtools import svg2paths
    paths, _ = svg2paths(os.path.join(HERE, "web", "img", "brand", "byk-interpretace.svg"))
    polys = []
    for p in paths:
        pts = [p.point(t) for t in np.linspace(0, 1, 300)]
        polys.append(SPoly([(z.real, -z.imag) for z in pts]).buffer(0))
    g = unary_union(polys); minx, miny, maxx, maxy = g.bounds; s = width_mm / (maxx - minx)
    from shapely import affinity
    return affinity.scale(affinity.translate(g, -minx, -miny), s, s, origin=(0, 0))
def emboss_y(poly2d, y_lo, y_hi, x_off, z_off, mirror=False):
    """polygony v rovině (x,z) vytáhnout podél y mezi y_lo a y_hi (skutečné souřadnice), posunout o (x_off, z_off);
    mirror zrcadlí text podle svislé osy (pro levobok), zrcadlení je kolem x=0 polygonu, pak se posune."""
    from shapely import affinity
    g = affinity.scale(poly2d, -1, 1, origin=(0, 0)) if mirror else poly2d
    geoms = list(g.geoms) if hasattr(g, "geoms") else [g]
    h = abs(y_hi - y_lo); meshes = []
    for pg in geoms:
        if pg.is_empty or pg.area < 0.05: continue
        m = extrude_polygon(pg, h)                                                   # lokální z = vytažení
        m.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, (1, 0, 0)))   # (x,Y,Z) -> (x,-Z,Y): vytažení do -y
        m.apply_translation((x_off, y_hi, z_off))                                    # y v [y_lo, y_hi]
        meshes.append(m)
    out = U(meshes) if len(meshes) > 1 else meshes[0]; out.fix_normals(); return out

# ---------------- pevná část zádě: deska a pouzdra čepů kormidel
rud_x = 10.5; rud_y = 5.5; rud_zk = 14.5
stab_plate = B(12.0, 17.0, 1.6, 10.0, 0, 20.3)
stab_boss = [Cz(6.0, 5.2, rud_x, s * rud_y, 17.1) for s in (-1, 1)]
solid = U([hull, stab_plate] + stab_boss)
deck_z = ZD * SCALE

def prism_xz(pts, y0, y1):
    """polygon v rovině (x, z) vytažený podél y"""
    m = extrude_polygon(SPoly(pts), y1 - y0)
    m.apply_transform(np.array([[1, 0, 0, 0], [0, 0, 1, y0], [0, 1, 0, 0], [0, 0, 0, 1]], float))
    if m.volume < 0: m.invert()
    return m
def frust(l0, w0, z0, l1, w1, z1, cx):
    pts = [(cx + sx * l0 / 2, sy * w0 / 2, z0) for sx in (-1, 1) for sy in (-1, 1)] + [(cx + sx * l1 / 2, sy * w1 / 2, z1) for sx in (-1, 1) for sy in (-1, 1)]
    return trimesh.convex.convex_hull(np.array(pts))

# ---------------- komora a sifon
cup = prism_xz([(ch_x0, -1), (ch_x1, -1), (ch_x1, CH_TOP_F), (ch_x0, CH_TOP_R)], -CH_W / 2, CH_W / 2)
notch = B(4.0, CH_W + 4.0, 2.4, ch_x0 - 1.5, 0, 0.0)
z_top = CH_TOP_F; z_exit = EXIT_Z
siphon = U([
    B(xb - ch_x1 + 0.1, SIPHON_D, PORT_H, (ch_x1 - 0.1 + xb) / 2, 0, CH_TOP_F - PORT_H / 2),      # okno u stropu komory
    Cz(SIPHON_D, CH_TOP_F - KEEL_T, xb, 0, (CH_TOP_F + KEEL_T) / 2),                          # sestupná větev
    B(xc - xb, SIPHON_D, WEIR_Z - KEEL_T, (xb + xc) / 2, 0, (WEIR_Z + KEEL_T) / 2),             # přepad: strop průchodu = WEIR_Z
    Cz(SIPHON_D, EXIT_Z + r - KEEL_T, xc, 0, (EXIT_Z + r + KEEL_T) / 2),                      # vzestupná větev
    Cy(SIPHON_D, 30.0, xc, 15.0, EXIT_Z),                                                    # výstup na pravobok
])
# ---------------- zátěž: matice na hraně ve štěrbinách kýlu
def nut_slot(x, n, stack=1):
    """šestihranná kapsa na matice naplocho (rohy podél osy lodi), otevřená zespodu; dvě lisovací žebírka na bočních ploškách, zářez na vypáčení"""
    zk_ = float(f_zk(x / SCALE)) * SCALE
    R = (n["e"] + 0.4) / 2; af_p = 2 * R * math.cos(math.radians(30)); dp = stack * n["t"] + 0.4; rib = 0.4
    hexp = trimesh.creation.cylinder(radius=R, height=dp + 1.0, sections=6); hexp.apply_translation((x, 0, zk_ - 1.0 + (dp + 1.0) / 2))
    pry = B(2.4, 4.0, 4.0, x - R - 0.8, 0, zk_ - 1.0 + 2.0)                                       # zářez na vypáčení šroubovákem
    ribs = [B(0.8, rib + 0.2, dp - 1.0, x + dx, s * (af_p / 2 - rib / 2 + 0.1), zk_ + 1.0 + (dp - 1.0) / 2) for dx in (-1.6, 1.6) for s in (-1, 1)]
    return [hexp, pry], ribs, np.array([x, 0.0, zk_ + dp - 0.1 - n["t"] / 2])                      # poloha první (nejhlouběji zatlačené) matice
slot_cav = []; slot_ribs = []; NUT_POS = {}
for x in M8_X:
    c_, r_, p_ = nut_slot(x, M8); slot_cav += c_; slot_ribs += r_; NUT_POS[("M8", x)] = p_
for x in M4_X:
    c_, r_, p_ = nut_slot(x, M4, stack=2); slot_cav += c_; slot_ribs += r_; NUT_POS[("M4", x)] = p_
# ---------------- šachta paliva pod věží: rovná, věž dosedá na palubu (hladké plochy z podložky), západky ji přitlačí
shaft = B(FUEL_L + 0.4, FUEL_W + 0.4, deck_z + 2 - (CH_TOP_R - 1.0), FUEL_X, 0, (deck_z + 2 + CH_TOP_R - 1.0) / 2)
BLADE_TOP = deck_z
HOOK = 0.6; ARM_L, ARM_T, ARM_W = 8.5, 1.0, FUEL_W
arm_bot = 21.5; arm_top = arm_bot + ARM_L; HOOK_Z0 = arm_bot + 0.2
HOOK_ZR = HOOK_Z0 + HOOK / math.tan(math.radians(30)); HOOK_ZF = HOOK_ZR + 0.4; HOOK_ZT = HOOK_ZF + HOOK
PRELOAD = 0.2                                                 # jazýček zůstane po zacvaknutí prohnutý o 0,2 mm a šikmou plochou zobáčku tlačí věž k palubě
DIMPLE_TOP = HOOK_ZT - (PRELOAD + 0.2)                        # horní hrana důlku: zobáček tam vyčnívá o PRELOAD víc, než je vůle šachty (0,2)
assert arm_top <= BLADE_TOP - 1.0, (arm_top, BLADE_TOP)
snap_dimples = [B(1.4, FUEL_W + 1.0, DIMPLE_TOP - (HOOK_Z0 - 0.15), FUEL_X + s * (FUEL_L / 2 + 0.2 + 0.5), 0, (HOOK_Z0 - 0.15 + DIMPLE_TOP) / 2) for s in (-1, 1)]
gun_x = (L_REAL - 26.0) * SCALE; flak_x = (L_REAL - 40.0) * SCALE
gun_hole = Cz(4.3, 0.8, gun_x, 0, deck_z - 0.2); flak_hole = Cz(3.3, 0.8, flak_x, 0, deck_z - 0.2)
# ---------------- hloubková kormidla: průchozí otvor pro hřídel, dosedací plochy, důlky aretace, rysky
plane_bow_x, plane_bow_z = 60.5 * SCALE, 3.3 * SCALE
plane_stern_x, plane_stern_z = 4.0 * SCALE, 3.4 * SCALE
plane_holes = [Cy(PLANE_SHAFT_D + 0.3, 60, plane_bow_x, 0, plane_bow_z), Cy(PLANE_SHAFT_D + 0.3, 60, plane_stern_x, 0, plane_stern_z)]
dimples = []; spot_faces = []
for px, pz in ((plane_bow_x, plane_bow_z), (plane_stern_x, plane_stern_z)):
    hw0 = half_width(px / SCALE, pz)
    for s in (-1, 1): spot_faces.append(Cy(12.8, 20.0, px, s * (hw0 + 10.0), pz))
    for ang in (-20, -10, 0, 10, 20):
        a = math.radians(ang); rr = 4.6
        mat = half_width((px + rr * math.sin(a)) / SCALE, pz + rr * math.cos(a)) - hw0
        assert mat > 0.1, f"důlek {ang}° u x={px:.0f}: trup nesahá k dosedací ploše ({mat:.2f} mm)"
        for s in (-1, 1): dimples.append(Sph(1.3, px + rr * math.sin(a), s * hw0, pz + rr * math.cos(a)))
rudder_holes = [Cz(RUDDER_PIN_D + 0.15, 7.0, rud_x, s * rud_y, rud_zk + 3.0) for s in (-1, 1)]
PIN_Y, PIN_Z = 4.5, 12.0
pin_holes = [Cx(PIN_D + 0.15, PIN_L + 2.0, SPLIT_X, s * PIN_Y, PIN_Z) for s in (-1, 1)]
hatches = [Cz(4.0, 1.0, x, 0, deck_z) for x in (12.0 * SCALE, 17.5 * SCALE, 55.0 * SCALE)]
slots = []
for x in np.arange(7.0, 63.5, 1.4):
    y = float(f_wc(x)) * SCALE
    for s in (-1, 1): slots.append(B(4.0, 1.0, 1.3, x * SCALE, s * y, deck_z - 2.6))      # mělké (0,5 mm): za nimi je stěna vzduchové dutiny
ticks = []
for px, pz in ((plane_bow_x, plane_bow_z), (plane_stern_x, plane_stern_z)):
    for ang in (-20, -10, 0, 10, 20):
        a = math.radians(ang); rr = 6.4
        for s in (-1, 1):
            yw = half_width(px / SCALE, pz + rr * math.cos(a)) * s
            ticks.append(B(0.5, 3.0, 1.4 if ang == 0 else 0.9, px + rr * math.sin(a), yw, pz + rr * math.cos(a)))
open_cav = U([cup, notch, siphon, shaft] + snap_dimples + [gun_hole, flak_hole] + slot_cav + plane_holes + spot_faces + dimples + rudder_holes + pin_holes + hatches + slots + ticks)

# ---------------- reliéf WILHELM na bocích
name_h = 5.0; name_z = 18.0; name_x0 = 46.0
name_poly = text_polys("WILHELM", name_h); nb = name_poly.bounds
emb = []
letters = list(name_poly.geoms) if hasattr(name_poly, "geoms") else [name_poly]
for s in (1, -1):
    for lp in letters:
        lb = lp.bounds; xc_l = (lb[0] + lb[2]) / 2
        x_world = (name_x0 + nb[2] - xc_l) if s == 1 else (name_x0 - nb[0] + xc_l)
        yw = half_width(x_world / SCALE, name_z + name_h / 2)
        if s == 1: emb.append(emboss_y(lp, yw - 1.6, yw + 0.6, name_x0 + nb[2], name_z, mirror=True))
        else:      emb.append(emboss_y(lp, -(yw + 0.6), -(yw - 1.6), name_x0 - nb[0], name_z))
solid = U([solid] + emb)
A = U([D(solid, open_cav)] + slot_ribs)
assert A.is_watertight
V_A = A.volume / 1000.0; C_A = A.center_mass.copy()
# kontrola: štěrbiny zátěže a komora nikde neprorazí bok trupu (stěna aspoň 1,2 mm)
def _grown(c):                                               # kapsa zvětšená o 1,2 mm do stran, od 0,4 mm nad místním kýlem
    x0_, x1_ = c.bounds[0][0] - 1.2, c.bounds[1][0] + 1.2; zb0 = max(float(f_zk(x0_ / SCALE)), float(f_zk(x1_ / SCALE))) * SCALE + 0.4
    return B(x1_ - x0_, c.extents[1] + 2.4, c.bounds[1][2] - zb0, (x0_ + x1_) / 2, 0, (c.bounds[1][2] + zb0) / 2)
_grow = U([_grown(c) for c in slot_cav[::2]] + [B(CH_L + 3.0, CH_W + 3.0, CH_TOP_R, cup_cx, 0, CH_TOP_R / 2 + 0.5)])
_out = D(_grow, hull); _v = 0.0 if _out.is_empty else max(q.volume for q in _out.split(only_watertight=False))
if _v >= 0.6: print("mimo trup:", [(q.bounds.round(1).tolist(), round(q.volume, 2)) for q in _out.split(only_watertight=False)])
assert _v < 0.6, f"štěrbina nebo komora je blíž než 1,2 mm k povrchu trupu ({_v:.1f} mm3 mimo)"

# ---------------- vzduchová dutina: tvar trupu zmenšený o stěnu, od roviny z_lo nahoru, s podélnou přepážkou
from shapely.geometry import LineString
def inner_ring(x_m, n=96):
    poly = SPoly(section(x_m, flat_fn(x_m))).buffer(0); inner = poly.buffer(-WALL)
    if inner.is_empty: return None
    wc_ = float(f_wc(x_m)) * SCALE                                    # rohy pod palubou vynechat: zvenku jsou tam štěrbiny v boku nástavby
    inner = inner.difference(SPoly([(wc_ - 2.4, 28.6), (40, 28.6), (40, 40), (wc_ - 2.4, 40)])).difference(SPoly([(-wc_ + 2.4, 28.6), (-40, 28.6), (-40, 40), (-wc_ + 2.4, 40)]))
    if inner.is_empty: return None
    if inner.geom_type != "Polygon": inner = max(inner.geoms, key=lambda g_: g_.area)
    zk_ = float(f_zk(x_m)) * SCALE; cz_ = zk_ + 0.5 * (deck_z - zk_)
    if inner.area < 8.0 or not inner.contains(Point(0, cz_)): return None
    ring = []
    for a in np.linspace(0, 2 * np.pi, n, endpoint=False):
        hit = LineString([(0, cz_), (80 * math.cos(a), cz_ + 80 * math.sin(a))]).intersection(inner.exterior)
        pts = [hit] if hit.geom_type == "Point" else list(getattr(hit, "geoms", []))
        d = min(math.hypot(p.x, p.y - cz_) for p in pts if p.geom_type == "Point")
        ring.append((d * math.cos(a), cz_ + d * math.sin(a)))
    ring = np.array(ring); ratio = SPoly(ring).area / inner.area
    return ring, cz_, ratio
rings = []
for x in stations:
    if not (2.4 <= x <= 64.2): continue
    rr_ = inner_ring(x)
    if rr_ is not None: rings.append((x * SCALE,) + rr_)
print(f"vnitřní trup: {len(rings)} řezů, x {rings[0][0]:.0f}–{rings[-1][0]:.0f} mm, nejhorší shoda paprskového řezu {min(q[3] for q in rings)*100:.1f} %")
def loft_rings(rs):
    M_ = len(rs[0][1]); V = []; F = []
    for x, ring, cz_, _ in rs: V += [(x, p[0], p[1]) for p in ring]
    for i in range(len(rs) - 1):
        a, b = i * M_, (i + 1) * M_
        for j in range(M_):
            j2 = (j + 1) % M_; F += [[a + j, b + j, b + j2], [a + j, b + j2, a + j2]]
    c0 = len(V); V.append((rs[0][0], 0, rs[0][2])); F += [[c0, (j + 1) % M_, j] for j in range(M_)]
    base = (len(rs) - 1) * M_; c1 = len(V); V.append((rs[-1][0], 0, rs[-1][2])); F += [[c1, base + j, base + (j + 1) % M_] for j in range(M_)]
    m = trimesh.Trimesh(np.array(V), np.array(F), process=True); m.fix_normals()
    if m.volume < 0: m.invert()
    return m
inner_hull = loft_rings(rings); assert inner_hull.is_volume, (inner_hull.is_watertight, inner_hull.is_winding_consistent, inner_hull.volume)
# plný materiál uvnitř obrysu dutiny. Trup se tiskne palubou dolů, proto každý blok buď stojí na palubě (sloup), nebo jde od boku k boku (most):
_gx0, _gx1 = ch_x0 - 2.0, xc + r + 2.0
keep = U([
    B(_gx1 - _gx0, CH_W + 4.0, 2 * (CH_TOP_F + 2.0), (_gx0 + _gx1) / 2, 0, 0),                             # komora a sifon + 2 mm stěna; po stranách zůstávají kapsy dutiny
    B(_gx1 - _gx0, 60, 1.6, (_gx0 + _gx1) / 2, 0, CH_TOP_F + 2.0 - 0.8),                                  # deska nad komorou od boku k boku (most), nad ní pokračuje dutina v nástavbě
    B(SIPHON_D + 4.0, 60, 2 * (EXIT_Z + r + 2.0), xc, 0, 0),                                              # vzestupná větev a výstup: od boku k boku
    B(FUEL_L + 0.4 + 2 * 1.4 + 4.0, FUEL_W + 0.4 + 4.0 + 1.0, 80, FUEL_X, 0, 20),                         # šachta paliva se západkami (sloup od paluby)
    B(PIN_L + 2.0 + 4.0, 60, 80, SPLIT_X, 0, 20),                                                         # dělicí rovina a kolíky
    Cz(8.0, 40, gun_x, 0, 20), Cz(7.0, 40, flak_x, 0, 20),
    B(17.0, 60, 80, plane_bow_x, 0, 20), B(17.0, 60, 80, plane_stern_x, 0, 20),                               # hřídele kormidel, dosedací plochy a rysky
    B(14.0, 60, 2 * 24.0, rud_x, 0, 0),
] + [B(n_["e"] + 0.4 + 4.0, 60, 2 * (k_ * n_["t"] + 0.4 + 2.0), x_, 0, float(f_zk(x_ / SCALE)) * SCALE) for n_, xs_, k_ in ((M8, M8_X, 1), (M4, M4_X, 2)) for x_ in xs_])   # kapsy zátěže: strop od boku k boku
VOID_XA = 14.0
def make_void(z_lo, x_f):
    v = trimesh.boolean.intersection([inner_hull, B(x_f - VOID_XA, 80, 60, (x_f + VOID_XA) / 2, 0, z_lo + 30)])
    return D(v, keep)

# ---------------- věž s U 47 a býkem
tower_h = (9.5 - 6.5) * SCALE; tw_l, tw_w = 6.4 * SCALE, 1.9 * SCALE
tower = extrude_polygon(stadium(tw_l, tw_w), tower_h)
well = extrude_polygon(stadium(tw_l - 3.2, tw_w - 3.2), 3.0); well.apply_translation((0, 0, tower_h - 3.0 + 0.01)); tower = D(tower, well)
scopes = [Cz(2.2, 9.0, -4.0, 1.2, tower_h - 3.0 + 4.5), Cz(1.8, 7.0, -0.5, -1.4, tower_h - 3.0 + 3.5), Cz(1.6, 3.5, 9.0, 0, tower_h - 3.0 + 1.75)]
num_poly = text_polys("U 47", 3.6); nbb = num_poly.bounds
bull = bull_poly(6.5); bb = bull.bounds
tw_emb = []
for s in (1, -1):
    yw = tw_w / 2
    if s == 1:
        tw_emb.append(emboss_y(num_poly, yw - 0.6, yw + 0.4, -9.0 + nbb[2], 4.0, mirror=True))
        tw_emb.append(emboss_y(bull, yw - 0.6, yw + 0.4, 1.5 + (bb[2] - bb[0]), 3.5, mirror=True))
    else:
        tw_emb.append(emboss_y(num_poly, -(yw + 0.4), -(yw - 0.6), -9.0 - nbb[0], 4.0))
        tw_emb.append(emboss_y(bull, -(yw + 0.4), -(yw - 0.6), 1.5, 3.5))
tower = U([tower] + scopes + tw_emb)
tw_socket = B(4.2, FUEL_W + 0.4, 1.6, FUEL_X - tower_cx, 0, 0.8 - 0.01)       # drážka ve spodku věže pro výstupek kazety
assert tw_w - (FUEL_W + 0.4) >= 1.6, "věž je nad šachtou příliš úzká"
tower = D(tower, tw_socket)

# ---------------- palivová kazeta: čepel v šachtě, dva pružné jazýčky, dole košíček otevřený z boku (souřadnice x kolem 0, z skutečné)
BASKET_Z0, BASK_TOP = 2.9, 14.0
blade = B(FUEL_L, FUEL_W, BLADE_TOP - BASKET_Z0, 0, 0, (BLADE_TOP + BASKET_Z0) / 2)
bumps = []; arm_cuts = []
for s_ in (-1, 1):
    xo = s_ * (FUEL_L / 2)
    arm_cuts.append(B(0.8, FUEL_W + 2, ARM_L + 0.8, xo - s_ * (ARM_T + 0.4), 0, arm_bot - 0.8 + (ARM_L + 0.8) / 2))
    arm_cuts.append(B(ARM_T + 2.0, FUEL_W + 2, 0.8, xo - s_ * (ARM_T + 0.8) / 2 + s_ * 0.6, 0, arm_bot - 0.4))
    pts = []
    for yy in (-ARM_W / 2, ARM_W / 2):
        pts += [(xo - s_ * 0.4, yy, HOOK_Z0), (xo, yy, HOOK_Z0), (xo + s_ * HOOK, yy, HOOK_ZR), (xo + s_ * HOOK, yy, HOOK_ZF), (xo, yy, HOOK_ZT), (xo - s_ * 0.4, yy, HOOK_ZT)]
    bumps.append(trimesh.convex.convex_hull(np.array(pts)))
cavity = B(FUEL_L - 1.6, FUEL_W - 1.6, BASK_TOP - 0.4 - (BASKET_Z0 + 0.8), 0, 0, (BASK_TOP - 0.4 + BASKET_Z0 + 0.8) / 2)
window = B(FUEL_L - 3.6, 2.0, 13.2 - 6.5, 0, FUEL_W / 2, (13.2 + 6.5) / 2)                         # otevřený bok (+y): tudy se vkládá tableta
slots_b = [B(1.2, 2.0, 13.0 - 4.5, dx, -FUEL_W / 2, (13.0 + 4.5) / 2) for dx in (-4.0, 0.0, 4.0)]     # protější stěna: štěrbiny pro vodu
holes_b = [Cz(1.2, 3, dx, 0, BASKET_Z0 + 0.4) for dx in (-3.5, 0, 3.5)]
tongue = B(3.8, FUEL_W, 1.4, 0, 0, BLADE_TOP + 0.7)
cartridge = D(U([D(blade, U(arm_cuts)), tongue] + bumps), U([cavity, window] + slots_b + holes_b))
BASKET_CM3 = (FUEL_L - 1.6) * (FUEL_W - 1.6) * (BASK_TOP - 0.4 - BASKET_Z0 - 0.8) / 1000
_F = lambda w, t, L_, eps, E=3500: w * t * t * E * eps / (6 * L_)
mu = 0.3; _y = HOOK - 0.2; _eps = 1.5 * ARM_T * _y / ARM_L ** 2; _f = _F(ARM_W, ARM_T, ARM_L, _eps)
_fp = _F(ARM_W, ARM_T, ARM_L, 1.5 * ARM_T * PRELOAD / ARM_L ** 2)
print(f"západka: záběr {_y:.1f} mm, ohyb {_eps*100:.2f} %, zasunutí {2*_f*(mu+math.tan(math.radians(30)))/(1-mu*math.tan(math.radians(30))):.1f} N, vytažení {2*_f*(mu+1)/(1-mu):.1f} N, přítlak věže k palubě {2*_fp*(1-mu)/(1+mu):.1f} N")
_root_cut = U([B(ARM_T + 0.2, FUEL_W + 2, 0.6, s_ * (FUEL_L / 2 - ARM_T / 2), 0, arm_top - 0.2) for s_ in (-1, 1)])
_parts = D(cartridge, _root_cut).split(only_watertight=False)
assert len(_parts) == 3, f"jazýček kazety není volný (těles po odříznutí kořene: {len(_parts)})"
assert cartridge.is_watertight and len(cartridge.split()) == 1
gun = U([Cz(4.0, 2.5, 0, 0, 1.25), B(5.0, 3.0, 3.0, 0, 0, 4.0), Cx(1.6, 16.0, 9.0, 0, 4.6), B(2.6, 2.0, 1.6, -3.0, 0, 4.6)])
gun = trimesh.boolean.intersection([gun, B(60, 60, 5.3, 0, 0, 5.3 / 2)])
flak = U([Cz(3.0, 1.5, 0, 0, 0.75), Cz(1.4, 4.4, 0, 0, 3.7), Cx(1.0, 7.0, 3.0, 0, 5.4)])
flak = trimesh.boolean.intersection([flak, B(60, 60, 5.75, 0, 0, 5.75 / 2)])
pins = U([Cx(PIN_D - 0.1, PIN_L, 0, s * 6, PIN_D / 2) for s in (-1, 1)])

# ---------------- hloubková kormidla ze dvou dílů: A = hřídel + list, B = list s nábojem, který se nasune na konec hřídele a přilepí
HUB_B = 3.0
def plane_parts(x_at, z_at):
    hw = half_width(x_at / SCALE, z_at)
    def side(s, hub_t):
        return [Cy(6.0 if hub_t < 2 else 7.0, hub_t, 0, s * (hw + hub_t / 2), 0), B(1.4, 1.6, 3.6, 0, s * (hw + 0.8), 2.6 + 1.8), Sph(1.2, 0, s * (hw + 0.3), 4.6),
                B(PLANE_CHORD, PLANE_SPAN / 2 - 2.0, PLANE_T, 0, s * (hw + 1.6 + (PLANE_SPAN / 2 - 2.0) / 2), 0)]
    y0, y1 = -hw - 1.6, hw + HUB_B
    pa = U([Cy(PLANE_SHAFT_D - 0.15, y1 - y0, 0, (y0 + y1) / 2, 0)] + side(-1, 1.6))
    hole = trimesh.boolean.intersection([Cy(PLANE_SHAFT_D, HUB_B + 2.0, 0, hw + HUB_B / 2, 0), B(20, 20, 20, 0, hw, 10 - PLANE_T / 2 - 0.05)])   # otvor tvaru D
    pb = D(U(side(1, HUB_B)), hole)
    flat = lambda m: trimesh.boolean.intersection([m, B(300, 300, 300, 0, 0, 150 - PLANE_T / 2)])      # spodek srovnaný: díl leží na podložce
    return flat(pa), flat(pb)
plane_bow_a, plane_bow_b = plane_parts(plane_bow_x, plane_bow_z); plane_stern_a, plane_stern_b = plane_parts(plane_stern_x, plane_stern_z)
for m_ in (plane_bow_a, plane_bow_b, plane_stern_a, plane_stern_b): assert m_.is_watertight and len(m_.split()) == 1
def rudder_part():
    return U([B(2.0 * SCALE, 1.6, 2.6 * SCALE, 0, 0, -1.3 * SCALE), Cz(RUDDER_PIN_D - 0.05, 6.0, 0, 0, 3.0)])
rudder = trimesh.boolean.intersection([rudder_part(), B(300, 300, 300, 0, 150 - 0.8, 0)])
TUBE_H = 0.8 + 4 * M4["t"] + 0.4                           # trubička na 4 matice M4 naplocho
lids = {}
for name, wfit in (("lid_tight", CH_W + 0.1), ("lid_normal", CH_W - 0.05), ("lid_loose", CH_W - 0.25)):
    lid = B(CH_L - 0.4, wfit, PLUG_T, 0, 0, PLUG_T / 2)
    holes = [Cz(2.0, PLUG_T + 2, x, y, PLUG_T / 2) for x in np.linspace(-CH_L / 2 + 4, CH_L / 2 - 4, 5) for y in (-3.0, 3.0)]
    lribs = [B(5, 0.6, PLUG_T - 0.8, x, s * (wfit / 2), PLUG_T / 2) for x in (-9, 0, 9) for s in (-1, 1)]
    tube = D(Cz(M4["e"] + 0.3 + 2.0, TUBE_H, LID_TUBE_X, 0, PLUG_T + TUBE_H / 2 - 0.01), U([Cz(M4["e"] + 0.3, TUBE_H, LID_TUBE_X, 0, PLUG_T + 0.8 + TUBE_H / 2),
                                                                                         B(1.6, 20, TUBE_H, LID_TUBE_X, 0, PLUG_T + 2.0 + TUBE_H / 2)]))   # dno 0,8 mm nad víčkem, dvě svislé štěrbiny pro vodu
    lids[name] = D(U([lid, tube] + lribs), U(holes))
def nut_mesh(n, hole):
    m = trimesh.creation.cylinder(radius=n["e"] / 2, height=n["t"], sections=6); m = D(m, trimesh.creation.cylinder(radius=hole / 2, height=n["t"] + 2, sections=24))
    return m                                                    # naplocho: osa svisle, rohy podél x

# ---------------- sestava v poloze v lodi a bilance hmotnosti
tower_w = tower.copy(); tower_w.apply_translation((tower_cx, 0, deck_z))
cart_w = cartridge.copy(); cart_w.apply_translation((FUEL_X, 0, 0))
gun_w = gun.copy(); gun_w.apply_translation((gun_x, 0, deck_z)); flak_w = flak.copy(); flak_w.apply_translation((flak_x, 0, deck_z))
lid_w = lids["lid_normal"].copy(); lid_w.apply_translation((cup_cx, 0, 0))
planes_w = []
for (pa, pb), (px, pz) in (((plane_bow_a, plane_bow_b), (plane_bow_x, plane_bow_z)), ((plane_stern_a, plane_stern_b), (plane_stern_x, plane_stern_z))):
    for m_ in (pa, pb): q = m_.copy(); q.apply_translation((px, 0, pz)); planes_w.append(q)
rud_w = [rudder.copy().apply_translation((rud_x, s * rud_y, rud_zk)) for s in (-1, 1)]
pins_w = [Cx(PIN_D - 0.1, PIN_L, SPLIT_X, s * PIN_Y, PIN_Z) for s in (-1, 1)]
PLASTIC = [tower_w, cart_w, gun_w, flak_w, lid_w] + planes_w + rud_w + pins_w
LID_NUT_POS = [np.array([cup_cx + LID_TUBE_X, 0.0, PLUG_T + 0.8 + M4["t"] * (k_ + 0.5)]) for k_ in range(4)]
def _s4(x, k_): return ("kapsa x=%g" % x, NUT_POS[("M4", x)] - np.array([0, 0, M4["t"] * k_]), "flat")          # k_=0: první matice na dně kapsy
def _l4(k_): return ("víčko", LID_NUT_POS[k_], "flat")
# pořadí, v jakém se přidávají dolaďovací matice M4: střídavě za komoru a na víčko (sklon se skoro nemění), nakonec na záď a na příď
TRIM_SEQ = [_s4(M4_X[2], 0), _l4(0), _s4(M4_X[2], 1), _l4(1), _s4(M4_X[1], 0), _l4(2), _s4(M4_X[1], 1), _l4(3), _s4(M4_X[0], 0), _s4(M4_X[4], 0), _s4(M4_X[0], 1), _s4(M4_X[4], 1), _s4(M4_X[3], 0), _s4(M4_X[3], 1)]
def nuts(fill8, ntrim=NOMINAL_TRIM): return [(M8["g"], NUT_POS[("M8", x)]) for x in fill8] + [(M4["g"], p_) for _, p_, _ in TRIM_SEQ[:ntrim]]
def totals(void, fill8=M8_FILL, ntrim=NOMINAL_TRIM, gas=None):
    """hmotnost [g], těžiště, výtlak [cm3], střed vztlaku celé ponorky pod vodou; gas = (objem cm3, těžiště plynu)"""
    Vv = void.volume / 1000; Cv = void.center_mass
    mass = [((V_A - Vv) * PLA, (V_A * C_A - Vv * Cv) / (V_A - Vv))] + [(m_.volume / 1000 * PLA, m_.center_mass) for m_ in PLASTIC] + nuts(fill8, ntrim)
    disp = [(V_A, C_A)] + [(m_.volume / 1000, m_.center_mass) for m_ in PLASTIC] + [(g_ / STEEL, p_) for g_, p_ in nuts(fill8, ntrim)]
    if gas: disp.append(gas)
    Mt = sum(a_ for a_, _ in mass); cg = sum(a_ * np.array(c_) for a_, c_ in mass) / Mt
    Vt = sum(a_ for a_, _ in disp); cb_ = sum(a_ * np.array(c_) for a_, c_ in disp) / Vt
    return Mt, cg, Vt, cb_
# plyn v komoře těsně před vypuštěním (rovná poloha): komora nad přepadem bez kazety + sestupná větev + okno
_gz = B(400, 60, 60, 166, 0, WEIR_Z + 30)
gas_reg = D(trimesh.boolean.intersection([U([cup, siphon]), B(xb + r - ch_x0 + 1, 60, 60, (ch_x0 + xb + r) / 2, 0, WEIR_Z + 30)]), U([cart_w, lid_w]))
GAS_CM3 = gas_reg.volume / 1000; GAS_C = gas_reg.center_mass.copy()
print(f"plyn před vypuštěním {GAS_CM3:.2f} cm3, těžiště plynu x = {GAS_C[0]:.1f} mm; košíček {BASKET_CM3:.2f} cm3")
# řešení: rovina z_lo (tíha ve vodě) a přední konec dutiny x_f (podélné vyvážení). Cíl: prázdná mírně přídí dolů, plná stejně přídí nahoru.
def resid(z_lo, x_f):
    v = make_void(z_lo, x_f); Mt, cg, Vt, cb_ = totals(v)
    dx_target = 0.5 * GAS_CM3 * (GAS_C[0] - cb_[0]) / Mt
    return np.array([Mt - Vt - TARGET_W_G, (cg[0] - cb_[0]) - dx_target]), v
X_F = 318.0                                                  # dutina až do přídě; podélné vyvážení se dělá polohou matic
za, zb_ = 9.5, 18.0; (fa, _), (fb, _) = resid(za, X_F), resid(zb_, X_F)
if fa[0] > 0: print(f"POZOR: i s největší dutinou je ponorka o {fa[0]:.2f} g těžší než cíl"); p = np.array([za, X_F])
else:
    for it in range(12):
        zn = float(np.clip(za + (0 - fa[0]) * (zb_ - za) / (fb[0] - fa[0]), 9.5, 27.0)); fn, _ = resid(zn, X_F)
        za, fa, zb_, fb = zb_, fb, zn, fn
        if abs(fn[0]) < 0.02: break
    p = np.array([zb_, X_F])
f0, void = resid(*p); print(f"  z_lo {p[0]:.2f}: odchylka tíhy {f0[0]:+.3f} g, těžiště je o {f0[1]:+.3f} mm před cílem")
z_void, x_void = float(p[0]), float(p[1])
hull_final = D(A, void); assert hull_final.is_watertight
assert abs(hull_final.volume / 1000 - (V_A - void.volume / 1000)) < 0.02, "dutina zasahuje mimo plný materiál"
import trimesh.proximity as _prox
_s = void.sample(6000); _d = np.abs(_prox.signed_distance(A, _s)); print(f"nejmenší stěna kolem dutiny {_d.min():.2f} mm v bodě {_s[_d.argmin()].round(1)}")
assert _d.min() > 1.25, "stěna kolem dutiny je tenčí než 1,25 mm"
Mt, cm, Vt, cb = totals(void); W_NET = Mt - Vt; BG = cb[2] - cm[2]
Mg, cmg, Vg, cbg = totals(void, gas=(GAS_CM3, GAS_C))
FS = 0.001 * CH_W * CH_L ** 3 / 12 / Mt                       # ztráta stability volnou hladinou plynu [mm]
pitch = lambda dxx, bg: math.degrees(math.atan2(dxx, bg))
print(f"dutina od z={z_void:.1f} mm po x={x_void:.0f} mm, {void.volume/1000:.1f} cm3; hmotnost {Mt:.1f} g, výtlak {Vt:.1f} cm3, tíha ve vodě {W_NET:+.2f} g (s plynem {Mg - Vg:+.2f} g)")
print(f"těžiště z={cm[2]:.2f} x={cm[0]:.2f}; střed vztlaku z={cb[2]:.2f} x={cb[0]:.2f}; BG={BG:.2f} mm (volná hladina −{FS:.2f}); sklon prázdná {pitch(cm[0]-cb[0], BG):+.1f}°, plná {pitch(cmg[0]-cbg[0], cbg[2]-cmg[2]-FS):+.1f}° (kladný = příď dolů)")

# ---------------- dělení trupu a export
big = 1000
stern = trimesh.boolean.intersection([hull_final, B(SPLIT_X + 1, big, big, (SPLIT_X + 1) / 2 - 1, 0, 0)])
bow = trimesh.boolean.intersection([hull_final, B(big, big, big, SPLIT_X + big / 2, 0, 0)])
print("díly vodotěsné:", stern.is_watertight, bow.is_watertight, "délky:", round(stern.extents[0], 1), round(bow.extents[0], 1))
for n_, m_ in (("záď", stern), ("příď", bow)):
    for q in m_.split(only_watertight=False): print("   těleso", n_, round(q.volume / 1000, 3), "cm3", q.bounds.round(1).tolist())
assert stern.is_watertight and bow.is_watertight
def export(mesh, name, flip=False, rot=None):
    m = mesh.copy()
    if flip: m.apply_transform(trimesh.transformations.rotation_matrix(np.pi, (1, 0, 0)))
    if rot is not None: m.apply_transform(rot)
    m.apply_translation(-m.bounds[0])
    m0 = m.copy()
    m.vertices = np.round(m.vertices, 4); m.merge_vertices(); m.update_faces(m.nondegenerate_faces()); m.update_faces(m.unique_faces()); m.remove_unreferenced_vertices()
    if not m.is_watertight: m = m0
    if not m.is_watertight: print("  POZOR: díl není uzavřené těleso:", name)
    m.export(os.path.join(OUT, name + ".stl")); return m
RX = lambda a: trimesh.transformations.rotation_matrix(a, (1, 0, 0))
# odlamovací žebra pod okraji stabilizační desky (jen pro tisk)
_zt = stern.bounds[1][2]; _zp = 20.3 + 0.8; _ry = 17.0 / 2 - 0.6
ribs = U([B(12.0, 1.2, _zt - (_zp + 0.4), 10.0, s_ * _ry, (_zt + _zp + 0.4) / 2) for s_ in (-1, 1)] + [B(12.0, 0.5, 0.42, 10.0, s_ * _ry, _zp + 0.2) for s_ in (-1, 1)])
assert trimesh.boolean.intersection([ribs, hull_final.copy()]).volume < 0.12, "žebro koliduje s trupem"
stern_print = U([stern, ribs])
for f in os.listdir(OUT):
    if f.endswith(".stl"): os.remove(os.path.join(OUT, f))
export(stern_print, "hull_stern_print", flip=True); export(bow, "hull_bow_print", flip=True)
export(tower, "tower"); export(cartridge, "fuel_cartridge", rot=RX(np.pi / 2)); export(gun, "gun_88", rot=RX(np.pi)); export(flak, "flak_20", rot=RX(np.pi)); export(pins, "pins")
export(plane_bow_a, "plane_bow_a"); export(plane_bow_b, "plane_bow_b"); export(plane_stern_a, "plane_stern_a"); export(plane_stern_b, "plane_stern_b")
export(rudder, "rudder", rot=RX(np.pi / 2))
for name, lid in lids.items(): export(lid, name)
# ---------------- výstavní díly: vložka do komory, zátka výstupu, stojánek
PLUG_H = CH_TOP_R - 0.5
pocket_c = FUEL_X - cup_cx
plug_body = B(CH_L - 0.3, CH_W - 0.15, PLUG_H, 0, 0, PLUG_H / 2)
pocket = B(FUEL_L + 1.2, CH_W + 2, PLUG_H, pocket_c, 0, BASKET_Z0 - 0.6 + PLUG_H / 2)
plug_ribs = [B(6, 0.6, PLUG_H - 2.0, x, s * ((CH_W - 0.15) / 2), PLUG_H / 2) for x in (5.0, 11.0) for s in (-1, 1)]
outer = D(U([plug_body] + plug_ribs), pocket)
inner = D(B(CH_L - 0.3 - 2.4, CH_W - 0.15 - 2.4, PLUG_H - 2.4, 0, 0, PLUG_H / 2), B(FUEL_L + 1.2 + 2.4, CH_W + 2, PLUG_H, pocket_c, 0, BASKET_Z0 - 0.6 - 1.2 + PLUG_H / 2))
display_plug = D(outer, inner); float_air = inner.volume / 1000
FLOAT_GAIN = outer.volume / 1000 - display_plug.volume / 1000 * PLA - lids["lid_normal"].volume / 1000 * (1 - PLA)
cap_body = Cy(SIPHON_D - 0.35, 5.0, 0, -2.5, 0)
cap_ribs = [trimesh.creation.cylinder(radius=0.35, height=4.0, sections=12) for _ in range(3)]
for k, rb in enumerate(cap_ribs):
    a = 2 * math.pi * k / 3; rb.apply_transform(RX(np.pi / 2)); rb.apply_translation(((SIPHON_D / 2 - 0.2) * math.cos(a), -2.6, (SIPHON_D / 2 - 0.2) * math.sin(a)))
exit_cap = D(U([cap_body] + cap_ribs), B(SIPHON_D, 1.2, 1.0, 0, -0.5, 0))
export(exit_cap, "exit_cap", rot=RX(np.pi / 2)); export(display_plug, "display_plug")
st_x0, st_x1 = 70.0, 270.0
stand = D(U([B(st_x1 - st_x0, 44.0, 4.0, (st_x0 + st_x1) / 2, 0, -8.0)] + [B(10.0, 30.0, 16.0, x, 0, -2.0) for x in (100.0, 240.0)]), hull)
lab = text_polys("U-47  WILHELM", 7.0); lb = lab.bounds
lab_m = U([extrude_polygon(g_, 0.8) for g_ in (list(lab.geoms) if hasattr(lab, "geoms") else [lab]) if g_.area > 0.05])
lab_m.apply_translation(((st_x0 + st_x1) / 2 - (lb[0] + lb[2]) / 2, -18.0 - lb[1], -6.0))
stand = U([stand, lab_m]); stand.export(os.path.join(WEB, "stand.stl")); export(stand, "stand")
# zkušební kostka: stejná komora a sifon (bez šachty); prášek se dá pod víčko
tb_x0 = ch_x0 - 6.0; tb_x1 = xc + r + 4.0; tb_h = EXIT_Z + r + 4.0
test_block = D(B(tb_x1 - tb_x0, CH_W + 8.0, tb_h, (tb_x0 + tb_x1) / 2, 0, tb_h / 2), U([cup, notch, siphon]))
export(test_block, "test_block_print", flip=True); test_block.export(os.path.join(WEB, "test_block.stl"))
print("výstavní vložka", [round(v, 1) for v in display_plug.extents], f"vztlak navíc {FLOAT_GAIN:+.1f} g; zkušební kostka", [round(v, 1) for v in test_block.extents], round(test_block.volume / 1000 * PLA, 1), "g")

# ---------------- modely pro web (sestava v poloze v lodi)
for f in ("ballast_caps.stl",):
    if os.path.exists(os.path.join(WEB, f)): os.remove(os.path.join(WEB, f))
for mesh, name in ((stern, "hull_stern"), (bow, "hull_bow"), (void, "void"), (cup, "cup"), (siphon, "siphon"), (tower_w, "tower"), (cart_w, "cartridge"), (gun_w, "gun"), (flak_w, "flak"), (lid_w, "lid")):
    mesh.export(os.path.join(WEB, name + ".stl"))
ec = exit_cap.copy(); ec.apply_translation((xc, half_width(xc / SCALE, z_exit) - 2.5, z_exit)); ec.export(os.path.join(WEB, "exit_cap.stl"))
U(planes_w).export(os.path.join(WEB, "planes.stl")); U(rud_w).export(os.path.join(WEB, "rudders.stl"))
_n8 = nut_mesh(M8, 6.6); _n4 = nut_mesh(M4, 3.3)
_n4f = _n4
U([_n8.copy().apply_translation(NUT_POS[("M8", x)]) for x in M8_FILL] + [(_n4 if o_ == "edge" else _n4f).copy().apply_translation(p_) for _, p_, o_ in TRIM_SEQ[:NOMINAL_TRIM]]).export(os.path.join(WEB, "ballast.stl"))

ASM = os.path.join(HERE, "sestava_v4"); os.makedirs(ASM, exist_ok=True)
for n_, m_ in (("trup", hull_final), ("vez", tower_w), ("kazeta", cart_w), ("delo", gun_w), ("kulomet", flak_w), ("vicko", lid_w), ("hk_prid_a", planes_w[0]), ("hk_prid_b", planes_w[1]), ("hk_zad_a", planes_w[2]), ("hk_zad_b", planes_w[3]),
               ("kormidlo_l", rud_w[0]), ("kormidlo_p", rud_w[1]), ("kolik_l", pins_w[0]), ("kolik_p", pins_w[1])): m_.export(os.path.join(ASM, n_ + ".stl"))
for i_, x in enumerate(M8_X): _n8.copy().apply_translation(NUT_POS[("M8", x)]).export(os.path.join(ASM, f"m8_{i_}.stl"))
for i_, (_, p_, o_) in enumerate(TRIM_SEQ): (_n4 if o_ == "edge" else _n4f).copy().apply_translation(p_).export(os.path.join(ASM, f"m4_{i_:02d}.stl"))

def sec_paths(m, origin, normal, ax):
    s = m.section(plane_origin=origin, plane_normal=normal); out = []
    for e in s.discrete:
        pts = e[:, ax]; keep_ = [pts[0]]
        for p_ in pts[1:]:
            if np.linalg.norm(p_ - keep_[-1]) > 0.2: keep_.append(p_)
        out.append([[round(float(p_[0]), 2), round(float(p_[1]), 2)] for p_ in keep_])
    return out
json.dump({"long": sec_paths(hull_final, (0, 0.6 + 0.7, 0), (0, 1, 0), [0, 2]), "cross_cup": sec_paths(hull_final, (cup_cx + 8, 0, 0), (1, 0, 0), [1, 2]),
           "cross_siphon": sec_paths(hull_final, (xb, 0, 0), (1, 0, 0), [1, 2])}, open(os.path.join(HERE, "v4_sections.json"), "w"))
g = lambda m: round(m.volume / 1000 * PLA, 2)
n8, n4 = len(M8_FILL), NOMINAL_TRIM
stats = {"version": 4, "scale": 200, "length_mm": round(L, 1), "beam_mm": round(2 * 3.1 * SCALE, 1), "height_mm": round(9.5 * SCALE, 1), "deck_mm": round(deck_z, 1), "waterline_mm": round(4.7 * SCALE, 1),
         "displaced_cm3": round(Vt, 1), "void_cm3": round(void.volume / 1000, 1), "void_z_lo_mm": round(z_void, 1), "void_x_end_mm": round(x_void, 1),
         "mass_g": {"hull_stern": g(stern), "hull_bow": g(bow), "tower": g(tower), "cartridge": g(cartridge), "gun": g(gun), "flak": g(flak), "lid": g(lids["lid_normal"]),
                    "planes": round(sum(g(m_) for m_ in (plane_bow_a, plane_bow_b, plane_stern_a, plane_stern_b)), 2), "rudders": round(2 * rudder.volume / 1000 * PLA, 2), "pins": g(pins),
                    "steel": round(n8 * M8["g"] + n4 * M4["g"], 1), "test_block": g(test_block)},
         "mass_total_g": round(Mt, 1), "plastic_g": round(Mt - n8 * M8["g"] - n4 * M4["g"], 1),
         "weight_in_water_g": round(W_NET, 2), "weight_in_water_gas_g": round(Mg - Vg, 2), "gas_capacity_cm3": round(GAS_CM3, 2), "gas_x_mm": round(float(GAS_C[0]), 1),
         "cg_z_mm": round(float(cm[2]), 2), "cb_z_mm": round(float(cb[2]), 2), "cg_x_mm": round(float(cm[0]), 2), "cb_x_mm": round(float(cb[0]), 2), "bg_mm": round(float(BG), 2), "free_surface_mm": round(FS, 2),
         "ballast": {"m8_slots_x": M8_X, "m4_slots_x": M4_X, "m8_fill_x": M8_FILL, "nominal_trim": NOMINAL_TRIM, "m8_g": M8["g"], "m4_g": M4["g"], "steel_density": STEEL,
                     "m8_pos": {str(x): [round(float(v), 2) for v in NUT_POS[("M8", x)]] for x in M8_X},
                     "trim_seq": [{"where": w_, "pos": [round(float(v), 2) for v in p_]} for w_, p_, _ in TRIM_SEQ]},
         "pla_density": PLA, "target_w_g": TARGET_W_G, "pins": {"x": SPLIT_X, "y": PIN_Y, "z": PIN_Z, "d": PIN_D - 0.1, "l": PIN_L}, "weir_z": WEIR_Z, "plug_t": PLUG_T,
         "cup": {"w": CH_W, "l": CH_L, "top": CH_TOP_F, "top_rear": CH_TOP_R, "cx": round(cup_cx, 1), "x0": round(ch_x0, 1)},
         "siphon": {"d": SIPHON_D, "weir_z": WEIR_Z, "bend_z": WEIR_Z, "exit_z": EXIT_Z, "top_z": CH_TOP_F, "xa": round(xa, 1), "xb": round(xb, 1), "xc": round(xc, 1), "keel_t": KEEL_T, "port_h": PORT_H},
         "split_x": SPLIT_X, "fuel": {"x": round(FUEL_X, 1), "l": FUEL_L, "w": FUEL_W, "basket_cm3": round(BASKET_CM3, 2)}, "float_gain_g": round(FLOAT_GAIN, 1), "tower_cx": round(tower_cx, 1), "gun_x": round(gun_x, 1), "flak_x": round(flak_x, 1),
         "planes": {"bow_x": plane_bow_x, "bow_z": plane_bow_z, "stern_x": plane_stern_x, "stern_z": plane_stern_z, "span": PLANE_SPAN, "chord": PLANE_CHORD, "detents_deg": [-20, -10, 0, 10, 20]},
         "rudders": {"x": rud_x, "y": rud_y, "z": round(rud_zk, 1), "pin_d": RUDDER_PIN_D},
         "parts_extents_mm": {k: [round(v, 1) for v in m.extents] for k, m in (("hull_stern", stern), ("hull_bow", bow), ("tower", tower), ("test_block", test_block), ("plane_bow", plane_bow_a), ("rudder", rudder))}}
json.dump(stats, open(os.path.join(HERE, "v4_stats.json"), "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats["mass_g"], ensure_ascii=False)); print("hotovo:", OUT)
