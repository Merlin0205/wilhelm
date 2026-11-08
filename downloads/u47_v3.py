# -*- coding: utf-8 -*-
"""
Wilhelm v3 – kluzák na prášek: ponorka U-47 (typ VIIB) 1:200, 332,5 mm, dělený trup pro Bambu Lab P1S.

Co je nového proti v2 (podle simulace a rešerše):
- sifon Ø 9 mm (Eötvös 11: bubliny stoupají samy, ohyb se po vyfouknutí vždy znovu zalije), komora 14 mm, ohyb 7 mm
- komora posunutá 25 mm před střed vztlaku: plyn zvedá příď při stoupání, šroubky vpředu sklání příď při klesání
  -> s natočenými hloubkovými kormidly klouže vpřed na obou fázích (podvodní kluzák), kormidla zatáčí
- hloubková kormidla otočná na hřídeli Ø 4 se západkami 0/±10/±20°, kormidla otočná na čepech Ø 3
- reliéfní nápis WILHELM na bocích trupu, U 47 a funící býk na věži
- dělicí rovina posunuta na 240 mm (komora celá v záďovém dílu)

Spuštění: python u47_v3.py  -> stl_v3/, web/models/, v3_stats.json, v3_sections.json
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
PLA = 1.24; TARGET_NET_G = -1.5
SPLIT_X = 240.0
CUP_W, CUP_L, CUP_TOP = 12.0, 60.0, 14.0
CUP_OFFSET = 25.0                 # mm před objemovým těžištěm (moment plynu -> příď nahoru při stoupání)
PLUG_T = 2.4
SIPHON_D = 9.0; BEND_Z = 7.0
WALL_SIDE, WALL_DECK = 2.2, 2.0
PIN_D, PIN_L = 4.0, 24.0
PLANE_SHAFT_D = 4.0; PLANE_SPAN = 34.0; PLANE_CHORD = 9.0; PLANE_T = 1.6
RUDDER_PIN_D = 3.0
FONT = "/System/Library/Fonts/Supplemental/DIN Condensed Bold.ttf"

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "stl_v3"); os.makedirs(OUT, exist_ok=True)
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
        zb = 0.75; wb = keel_flat_half if z - zk_ < zb else keel_flat_half * max(0.0, 1 - (z - zk_ - zb) / (1.4 * zb)); w = max(w, wb)
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
cb_x = float(hull0.center_mass[0]); cup_cx = cb_x + CUP_OFFSET
tower_cx = (L_REAL - 33.8) * SCALE
print(f"objemové těžiště x = {cb_x:.1f} mm, komora x = {cup_cx:.1f} mm, věž x = {tower_cx:.1f} mm")
flat_half_m = (CUP_W / 2 + 2.0) / SCALE
def flat_fn(x_m):
    d = abs(x_m * SCALE - cup_cx) - CUP_L / 2
    return flat_half_m if d <= 0 else flat_half_m * max(0.0, 1 - d / 14.0)
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

# ---------------- kormidla pevná část: čepy; pohyblivá kormidla zvlášť
rud_x = 10.5; rud_y = 5.5; rud_zk = 14.5                      # kormidla visí z pouzder pod stabilizační deskou
stab_plate = B(12.0, 17.0, 1.6, 10.0, 0, 20.3)                 # vodorovná deska za trupem (skutečná VIIB má kormidla na stabilizátoru)
stab_boss = [Cz(6.0, 5.2, rud_x, s * rud_y, 17.1) for s in (-1, 1)]
solid = U([hull, stab_plate] + stab_boss)

# ---------------- otevřené dutiny
cup = B(CUP_L, CUP_W, CUP_TOP + 1, cup_cx, 0, (CUP_TOP + 1) / 2 - 1)
notch = B(4.0, CUP_W + 4.0, 2.4, cup_cx + CUP_L / 2 + 1.5, 0, 0.0)
r = SIPHON_D / 2
xa = cup_cx - CUP_L / 2 + 4.0 + r; xb = xa - (SIPHON_D + 3.0); xc = xb - (SIPHON_D + 3.0)
z_top = CUP_TOP + 3.0 + r; z_exit = CUP_TOP + 1.0 + r
siphon = U([
    Cz(SIPHON_D, z_top - CUP_TOP + 1, xa, 0, (CUP_TOP - 1 + z_top) / 2),
    Cx(SIPHON_D, xa - xb + SIPHON_D, (xa + xb) / 2, 0, z_top),
    Cz(SIPHON_D, z_top - BEND_Z + SIPHON_D, xb, 0, (z_top + BEND_Z) / 2),
    Cx(SIPHON_D, xb - xc + SIPHON_D, (xb + xc) / 2, 0, BEND_Z),
    Cz(SIPHON_D, z_exit - BEND_Z + SIPHON_D, xc, 0, (z_exit + BEND_Z) / 2),
    Cy(SIPHON_D, 30.0, xc, 15.0, z_exit),
])
screw_x = [60.0, 90.0, SPLIT_X + 14.0, SPLIT_X + 40.0]          # dva vzadu, dva vpředu (trim: příď dolů)
screws = [Cz(3.3, 13.0, x, 0, 5.5) for x in screw_x]
deck_z = ZD * SCALE
tower_recess = extrude_polygon(stadium(6.4 * SCALE + 0.4, 1.9 * SCALE + 0.4), 0.8); tower_recess.apply_translation((tower_cx, 0, deck_z - 0.4))
gun_x = (L_REAL - 26.0) * SCALE; flak_x = (L_REAL - 40.0) * SCALE
gun_hole = Cz(2.3, 5.0, gun_x, 0, deck_z - 2.0); flak_hole = Cz(1.8, 4.0, flak_x, 0, deck_z - 1.5)
# otočná hloubková kormidla: průchozí otvor pro hřídel
plane_bow_x, plane_bow_z = 60.5 * SCALE, 3.3 * SCALE
plane_stern_x, plane_stern_z = 4.0 * SCALE, 3.4 * SCALE
plane_holes = [Cy(PLANE_SHAFT_D + 0.3, 60, plane_bow_x, 0, plane_bow_z), Cy(PLANE_SHAFT_D + 0.3, 60, plane_stern_x, 0, plane_stern_z)]
# západkové důlky okolo otvoru na boku trupu: 0, ±10, ±20° (důlky v trupu, výstupek na náboji)
dimples = []
for px, pz in ((plane_bow_x, plane_bow_z), (plane_stern_x, plane_stern_z)):
    for ang in (-20, -10, 0, 10, 20):
        a = math.radians(ang); rr = 4.6
        for s in (-1, 1):
            yw = half_width(px / SCALE, pz + rr * math.cos(a)) * s      # povrch ve výšce důlku
            dimples.append(Sph(1.3, px + rr * math.sin(a), yw, pz + rr * math.cos(a)))
# kormidla otočná: čepy Ø3 zespodu do zádě
rudder_holes = [Cz(RUDDER_PIN_D + 0.15, 7.0, rud_x, s * rud_y, rud_zk + 3.0) for s in (-1, 1)]   # z 14 až 21
pin_holes = [Cx(PIN_D + 0.15, PIN_L + 2.0, SPLIT_X, s * 4.5, 12.0) for s in (-1, 1)]
hatches = [Cz(4.0, 1.0, x, 0, deck_z) for x in (12.0 * SCALE, 17.5 * SCALE, 55.0 * SCALE)]
slots = []
for x in np.arange(7.0, 63.5, 1.4):
    y = float(f_wc(x)) * SCALE
    for s in (-1, 1): slots.append(B(4.0, 1.8, 1.3, x * SCALE, s * y, deck_z - 2.6))
# stupnice úhlu u hloubkových kormidel: rysky 0, ±10, ±20 (zahloubené 0,3 mm) – malé obdélníky nad důlky
ticks = []
for px, pz in ((plane_bow_x, plane_bow_z), (plane_stern_x, plane_stern_z)):
    for ang in (-20, -10, 0, 10, 20):
        a = math.radians(ang); rr = 6.4
        for s in (-1, 1):
            yw = half_width(px / SCALE, pz + rr * math.cos(a)) * s
            tk = B(0.5, 3.0, 1.4 if ang == 0 else 0.9, px + rr * math.sin(a), yw, pz + rr * math.cos(a)); ticks.append(tk)   # 3 mm na tloušťku: vždy otevřené ven
open_cav = U([cup, notch, siphon, tower_recess, gun_hole, flak_hole] + screws + plane_holes + dimples + rudder_holes + pin_holes + hatches + slots + ticks)

# ---------------- reliéf: WILHELM na bocích (sedlové nádrže), zrcadlově na levoboku
name_h = 5.0; name_z = 18.0; name_x0 = 46.0
name_poly = text_polys("WILHELM", name_h)
nb = name_poly.bounds; name_w = nb[2] - nb[0]
emb = []
letters = list(name_poly.geoms) if hasattr(name_poly, "geoms") else [name_poly]
for s in (1, -1):
    for lp in letters:
        lb = lp.bounds; xc_l = (lb[0] + lb[2]) / 2
        # pozice písmene na trupu: pravobok zrcadlově (čte se zleva doprava při pohledu z pravoboku)
        x_world = (name_x0 + nb[2] - xc_l) if s == 1 else (name_x0 - nb[0] + xc_l)
        yw = half_width(x_world / SCALE, name_z + name_h / 2)
        if s == 1: emb.append(emboss_y(lp, yw - 1.6, yw + 0.6, name_x0 + nb[2], name_z, mirror=True))
        else:      emb.append(emboss_y(lp, -(yw + 0.6), -(yw - 1.6), name_x0 - nb[0], name_z))
solid = U([solid] + emb)
A = D(solid, open_cav)
V_disp = A.volume / 1000.0

# ---------------- dutina plováku
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
        B(xa - xc + 2 * SIPHON_D + 6.0, 60, 2 * (z_top + r + 2.5), (xa + xc) / 2, 0, 0),   # sifon + 2,5 mm stěna nad ním
        B(6.0, 60, 60, SPLIT_X, 0, 20),
        Cz(8.0, 40, gun_x, 0, 20), Cz(7.0, 40, flak_x, 0, 20),
        B(12.0, 60, 2 * (plane_bow_z + PLANE_SHAFT_D / 2 + 2.5), plane_bow_x, 0, 0),      # hřídele kormidel
        B(12.0, 60, 2 * (plane_stern_z + PLANE_SHAFT_D / 2 + 2.5), plane_stern_x, 0, 0),
        B(14.0, 60, 2 * 24.0, rud_x, 0, 0),                                                  # záď: deska a čepy kormidel
    ])
    return D(v, keep)
def net_for(z_lo):
    v = make_void(z_lo); H = D(A, v); return H, v, H.volume / 1000.0 * PLA - V_disp

# ---------------- věž s U 47 a býkem
tower_h = (9.5 - 6.5) * SCALE; tw_l, tw_w = 6.4 * SCALE, 1.9 * SCALE
tower = extrude_polygon(stadium(tw_l, tw_w), tower_h)
well = extrude_polygon(stadium(tw_l - 3.2, tw_w - 3.2), 3.0); well.apply_translation((0, 0, tower_h - 3.0 + 0.01)); tower = D(tower, well)
scopes = [Cz(2.2, 9.0, -4.0, 1.2, tower_h - 3.0 + 4.5), Cz(1.8, 7.0, -0.5, -1.4, tower_h - 3.0 + 3.5), Cz(1.6, 3.5, 9.0, 0, tower_h - 3.0 + 1.75)]
num_poly = text_polys("U 47", 3.6); nbb = num_poly.bounds; num_w = nbb[2] - nbb[0]
bull = bull_poly(6.5); bb = bull.bounds
tw_emb = []
for s in (1, -1):
    yw = tw_w / 2
    # číslo vzadu (x záporné = záď), býk vpředu
    if s == 1:
        tw_emb.append(emboss_y(num_poly, yw - 0.6, yw + 0.4, -9.0 + nbb[2], 4.0, mirror=True))
        tw_emb.append(emboss_y(bull, yw - 0.6, yw + 0.4, 1.5 + (bb[2] - bb[0]), 3.5, mirror=True))
    else:
        tw_emb.append(emboss_y(num_poly, -(yw + 0.4), -(yw - 0.6), -9.0 - nbb[0], 4.0))
        tw_emb.append(emboss_y(bull, -(yw + 0.4), -(yw - 0.6), 1.5, 3.5))
tower = U([tower] + scopes + tw_emb)
gun = U([Cz(4.0, 2.5, 0, 0, 1.25), Cz(2.0, 3.0, 0, 0, -1.5 + 0.01), B(5.0, 3.0, 3.0, 0, 0, 4.0), Cx(1.6, 16.0, 9.0, 0, 4.6), B(2.6, 2.0, 1.6, -3.0, 0, 4.6)])
flak = U([Cz(3.0, 1.5, 0, 0, 0.75), Cz(1.5, 2.5, 0, 0, -1.25 + 0.01), Cz(1.4, 4.0, 0, 0, 3.5), Cx(1.0, 7.0, 3.0, 0, 5.4)])
pins = U([Cx(PIN_D - 0.1, PIN_L, 0, s * 6, PIN_D / 2) for s in (-1, 1)])

# ---------------- otočná hloubková kormidla (hřídel + 2 listy + náboje s výstupkem), stejný díl pro příď i záď
def plane_part(x_at, z_at):
    hw = half_width(x_at / SCALE, z_at)
    shaft_len = 2 * hw + 2 * 3.2                                   # hřídel prochází trupem, náboje venku
    shaft = Cy(PLANE_SHAFT_D - 0.15, shaft_len, 0, 0, 0)
    parts = [shaft]
    for s in (-1, 1):
        hub = Cy(6.0, 1.6, 0, s * (hw + 0.8), 0)                    # náboj přiléhá k boku trupu
        arm = B(1.4, 1.6, 3.6, 0, s * (hw + 0.8), 2.6 + 1.8)        # ukazatel úhlu, nese západkový výstupek
        bump = Sph(1.2, 0, s * (hw + 0.3), 4.6)                     # výstupek zapadá do důlků 0/±10/±20°
        parts.append(arm)
        blade = B(PLANE_CHORD, PLANE_SPAN / 2 - 2.0, PLANE_T, 0, s * (hw + 1.6 + (PLANE_SPAN / 2 - 2.0) / 2), 0)
        parts += [hub, bump, blade]
    return U(parts)
plane_bow = plane_part(plane_bow_x, plane_bow_z); plane_stern = plane_part(plane_stern_x, plane_stern_z)
# kormidlo: list 10 × 1,6 × 13 mm s čepem Ø 3 × 6 mm nahoře (do zádě), ukazatel úhlu na horní hraně
def rudder_part():
    blade = B(2.0 * SCALE, 1.6, 2.6 * SCALE, 0, 0, -1.3 * SCALE)
    pin = Cz(RUDDER_PIN_D - 0.05, 6.0, 0, 0, 3.0)
    return U([blade, pin])
rudder = rudder_part()

tower_net = tower.volume / 1000 * (PLA - 1.0); gun_net = (gun.volume + flak.volume) / 1000 * (PLA - 1.0)
plug_net = CUP_L * CUP_W * PLUG_T / 1000 * (PLA - 1.0) * 0.85
pins_net = pins.volume / 1000 * PLA
planes_net = (plane_bow.volume + plane_stern.volume) / 1000 * (PLA - 1.0) + 0.3   # hřídel v otvoru: otvor je dutina
rudders_net = 2 * rudder.volume / 1000 * (PLA - 1.0)
target = TARGET_NET_G - tower_net - gun_net - plug_net - pins_net - planes_net - rudders_net
z_a, z_b = 22.0, 27.0
H_a, v_a, n_a = net_for(z_a); H_b, v_b, n_b = net_for(z_b)
for _ in range(7):
    z_new = float(np.clip(z_a + (target - n_a) * (z_b - z_a) / (n_b - n_a), 14.0, deck_z - WALL_DECK - 2.0))
    H_n, v_n, n_n = net_for(z_new)
    z_a, n_a, H_a, v_a = z_b, n_b, H_b, v_b; z_b, n_b, H_b, v_b = z_new, n_n, H_n, v_n
    if abs(n_n - target) < 0.05: break
hull_final, void, z_void = H_b, v_b, z_b
net_total = n_b + tower_net + gun_net + plug_net + pins_net + planes_net + rudders_net
cup_cap = CUP_L * CUP_W * (CUP_TOP - BEND_Z) / 1000
cm, cb = hull_final.center_mass, A.center_mass
print(f"dutina od z={z_void:.1f} mm, {void.volume/1000:.1f} cm3; vytlačeno {V_disp:.1f} cm3; PLA trup {hull_final.volume/1000*PLA:.0f} g; bilance {net_total:+.2f} g; plyn {cup_cap:.1f} cm3")
print(f"těžiště z={cm[2]:.1f} x={cm[0]:.1f}, střed vztlaku z={cb[2]:.1f} x={cb[0]:.1f} mm")

big = 1000
stern = trimesh.boolean.intersection([hull_final, B(SPLIT_X + 1, big, big, (SPLIT_X + 1) / 2 - 1, 0, 0)])
bow = trimesh.boolean.intersection([hull_final, B(big, big, big, SPLIT_X + big / 2, 0, 0)])
print("díly vodotěsné:", stern.is_watertight, bow.is_watertight, "délky:", round(stern.extents[0], 1), round(bow.extents[0], 1))

def export(mesh, name, flip=False, rot=None):
    m = mesh.copy()
    if flip: m.apply_transform(trimesh.transformations.rotation_matrix(np.pi, (1, 0, 0)))
    if rot is not None: m.apply_transform(rot)
    m.apply_translation(-m.bounds[0]); m.export(os.path.join(OUT, name + ".stl")); return m
export(stern, "hull_stern_print", flip=True); export(bow, "hull_bow_print", flip=True)
export(tower, "tower"); export(gun, "gun_88"); export(flak, "flak_20"); export(pins, "pins")
export(plane_bow, "plane_bow"); export(plane_stern, "plane_stern")
export(rudder, "rudder", rot=trimesh.transformations.rotation_matrix(np.pi / 2, (1, 0, 0)))   # naležato
lids = {}
for name, wfit in (("lid_tight", CUP_W + 0.1), ("lid_normal", CUP_W - 0.05), ("lid_loose", CUP_W - 0.25)):
    lid = B(CUP_L - 0.4, wfit, PLUG_T, 0, 0, PLUG_T / 2)
    holes = [Cz(2.0, PLUG_T + 2, x, y, PLUG_T / 2) for x in np.linspace(-CUP_L / 2 + 5, CUP_L / 2 - 5, 10) for y in (-3.0, 3.0)]
    ribs = [B(8, 0.6, PLUG_T - 0.8, x, s * (wfit / 2), PLUG_T / 2) for x in (-20, 0, 20) for s in (-1, 1)]
    lid = D(U([lid] + ribs), U(holes)); lid.export(os.path.join(OUT, name + ".stl")); lids[name] = lid
# ---------------- výstavní verze: vložka zaplní komoru, stojánek s nápisem
plug_body = B(CUP_L - 0.3, CUP_W - 0.15, CUP_TOP - 0.3, 0, 0, (CUP_TOP - 0.3) / 2)
plug_ribs = [B(8, 0.6, CUP_TOP - 2.0, x, s * ((CUP_W - 0.15) / 2), (CUP_TOP - 0.3) / 2) for x in (-20, 0, 20) for s in (-1, 1)]
display_plug = U([plug_body] + plug_ribs)
export(display_plug, "display_plug")
st_x0, st_x1 = 70.0, 270.0
base = B(st_x1 - st_x0, 44.0, 4.0, (st_x0 + st_x1) / 2, 0, -8.0)
cradles = [B(10.0, 30.0, 16.0, x, 0, -2.0) for x in (100.0, 240.0)]
stand = D(U([base] + cradles), hull)
lab = text_polys("U-47  WILHELM", 7.0); lb = lab.bounds
lab_geoms = list(lab.geoms) if hasattr(lab, "geoms") else [lab]
lab_m = U([extrude_polygon(g, 0.8) for g in lab_geoms if g.area > 0.05])
lab_m.apply_translation(((st_x0 + st_x1) / 2 - (lb[0] + lb[2]) / 2, -18.0 - lb[1], -6.0))
stand = U([stand, lab_m])
stand.export(os.path.join(WEB, "stand.stl"))
export(stand, "stand")
print("výstavní vložka", [round(v, 1) for v in display_plug.extents], "stojánek", [round(v, 1) for v in stand.extents], "vodotěsné", display_plug.is_watertight, stand.is_watertight)

# zkušební kostka: stejná komora, sifon a zátka
tb_x0 = xc - r - 4.0; tb_x1 = cup_cx + CUP_L / 2 + 6.0
test_block = D(B(tb_x1 - tb_x0, CUP_W + 8.0, z_top + r + 4.0, (tb_x0 + tb_x1) / 2, 0, (z_top + r + 4.0) / 2), U([cup, notch, siphon]))
export(test_block, "test_block_print", flip=True)
test_block.export(os.path.join(WEB, "test_block.stl"))
print("zkušební kostka", [round(v, 1) for v in test_block.extents], "mm,", round(test_block.volume / 1000 * PLA, 1), "g, vodotěsná", test_block.is_watertight)

# ---------------- modely pro web (sestava)
for mesh, name in ((stern, "hull_stern"), (bow, "hull_bow"), (void, "void"), (cup, "cup"), (siphon, "siphon")): mesh.export(os.path.join(WEB, name + ".stl"))
tower.copy().apply_translation((tower_cx, 0, deck_z)).export(os.path.join(WEB, "tower.stl"))
gun.copy().apply_translation((gun_x, 0, deck_z)).export(os.path.join(WEB, "gun.stl"))
flak.copy().apply_translation((flak_x, 0, deck_z)).export(os.path.join(WEB, "flak.stl"))
lids["lid_normal"].copy().apply_translation((cup_cx, 0, 0)).export(os.path.join(WEB, "lid.stl"))
U([plane_bow.copy().apply_translation((plane_bow_x, 0, plane_bow_z)), plane_stern.copy().apply_translation((plane_stern_x, 0, plane_stern_z))]).export(os.path.join(WEB, "planes.stl"))
U([rudder.copy().apply_translation((rud_x, s * rud_y, rud_zk)) for s in (-1, 1)]).export(os.path.join(WEB, "rudders.stl"))

def sec_paths(m, origin, normal, ax):
    s = m.section(plane_origin=origin, plane_normal=normal); out = []
    for e in s.discrete:
        pts = e[:, ax]; keep = [pts[0]]
        for p in pts[1:]:
            if np.linalg.norm(p - keep[-1]) > 0.2: keep.append(p)
        out.append([[round(float(p[0]), 2), round(float(p[1]), 2)] for p in keep])
    return out
json.dump({"long": sec_paths(hull_final, (0, 0, 0), (0, 1, 0), [0, 2]), "cross_cup": sec_paths(hull_final, (cup_cx, 0, 0), (1, 0, 0), [1, 2]),
           "cross_siphon": sec_paths(hull_final, (xb, 0, 0), (1, 0, 0), [1, 2])}, open(os.path.join(HERE, "v3_sections.json"), "w"))
g = lambda m: round(m.volume / 1000 * PLA, 2)
stats = {"version": 3, "scale": 200, "length_mm": round(L, 1), "beam_mm": round(2 * 3.1 * SCALE, 1), "height_mm": round(9.5 * SCALE, 1), "deck_mm": round(deck_z, 1), "waterline_mm": round(4.7 * SCALE, 1),
         "displaced_cm3": round(V_disp, 1), "void_cm3": round(void.volume / 1000, 1), "void_z_lo_mm": round(z_void, 1),
         "mass_g": {"hull_stern": g(stern), "hull_bow": g(bow), "tower": g(tower), "gun": g(gun), "flak": g(flak), "lid": round(CUP_L * CUP_W * PLUG_T / 1000 * PLA * 0.85, 2), "plane_bow": g(plane_bow), "plane_stern": g(plane_stern), "rudders": round(2 * rudder.volume / 1000 * PLA, 2), "pins": g(pins), "test_block": g(test_block)},
         "net_g": round(net_total, 2), "gas_capacity_cm3": round(cup_cap, 1), "cg_z_mm": round(cm[2], 1), "cb_z_mm": round(cb[2], 1), "cg_x_mm": round(cm[0], 1), "cb_x_mm": round(cb[0], 1),
         "cup": {"w": CUP_W, "l": CUP_L, "top": CUP_TOP, "cx": round(cup_cx, 1), "offset_mm": CUP_OFFSET}, "siphon": {"d": SIPHON_D, "bend_z": BEND_Z, "exit_z": round(z_exit, 1), "top_z": round(z_top, 1), "xa": round(xa, 1), "xb": round(xb, 1), "xc": round(xc, 1)},
         "split_x": SPLIT_X, "tower_cx": round(tower_cx, 1), "gun_x": round(gun_x, 1), "flak_x": round(flak_x, 1), "screw_x": screw_x,
         "planes": {"bow_x": plane_bow_x, "bow_z": plane_bow_z, "stern_x": plane_stern_x, "stern_z": plane_stern_z, "span": PLANE_SPAN, "chord": PLANE_CHORD, "detents_deg": [-20, -10, 0, 10, 20]},
         "rudders": {"x": rud_x, "y": rud_y, "z": round(rud_zk, 1), "pin_d": RUDDER_PIN_D},
         "parts_extents_mm": {k: [round(v, 1) for v in m.extents] for k, m in (("hull_stern", stern), ("hull_bow", bow), ("tower", tower), ("test_block", test_block), ("plane_bow", plane_bow), ("rudder", rudder))}}
stats["mass_g"]["total"] = round(sum(v for k, v in stats["mass_g"].items() if k != "test_block"), 1)
json.dump(stats, open(os.path.join(HERE, "v3_stats.json"), "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats["mass_g"], ensure_ascii=False)); print("hotovo:", OUT)
