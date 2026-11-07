# -*- coding: utf-8 -*-
"""
Fyzikální simulátor hračky Wilhelm: komora na prášek + plynový sifon + vztlak + sklon + kluz.
Soustředěný (lumped) model, krok 1 ms. Všechny jednotky SI, geometrie se zadává v mm.

Co model obsahuje
- kinetika prášku do pečiva: rychlá a pomalá složka uvolňování CO2, reakce jen když je prášek pod vodou
- hydrostatika komory: plyn u stropu, hladina v komoře, tlak plynu = tlak vody na hladině
- kapilarita: plyn vstoupí do kanálku / otvorů zátky jen po překonání 2σ·cosθ/r (ustupující úhel vody)
- sifon: spuštění, výtok (setrvačnost vodního sloupce + odpor), ukončení a znovuzalití ohybu
  (kritérium Eötvösova čísla pro samovolný výstup bublin, nebo tlakový rozdíl vs. hystereze)
- únik dírkami zátky, když sifon nespustí
- svislý pohyb: vztlak, tíha, odpor (plošný), přidaná hmotnost; hladina a dno nádrže
- sklon: moment plynu před/za těžištěm vs. vratný moment (těžiště pod středem vztlaku), kvazistatický
- kluz: hloubková kormidla s nastavitelným úhlem, vztlak a odpor desky při malém Re, odpor trupu; poloha x
- kormidlo: odhad poloměru kruhu

Spuštění: python wilhelm_sim.py [--mc]   (vykreslí grafy do sim/out/ a vypíše souhrn)
"""
import math, json, os, sys, argparse
import numpy as np

RHO = 1000.0; G = 9.81; SIGMA = 0.072; RHO_GAS = 1.8; NU = 1.0e-6
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ návrh (mm, g)
def design_v2(stats="v2_stats.json"):
    st = json.load(open(os.path.join(HERE, "..", stats)))
    return dict(
        name="v2", L_hull=st["length_mm"], B_hull=st["beam_mm"], H_deck=st["deck_mm"],
        mass_print_g=st["mass_g"]["total"], disp_cm3=st["displaced_cm3"] + (st["mass_g"]["total"] - st["mass_g"]["hull_stern"] - st["mass_g"]["hull_bow"]) / 1.24, cg_z=st["cg_z_mm"], cb_z=st["cb_z_mm"],
        cup_w=st["cup"]["w"], cup_l=st["cup"]["l"], cup_top=st["cup"]["top"], plug_top=2.4, hole_d=2.0, n_holes=20,
        siphon_d=st["siphon"]["d"], bend_z=st["siphon"]["bend_z"], exit_z=st["siphon"]["exit_z"], inlet_z=st["cup"]["top"], disp_extra=0.0,
        cup_x_offset=st["cup"].get("offset_mm", 0.0),          # mm před středem vztlaku (kladné = k přídi)
        plane_area_mm2=(4 * (st["planes"]["span"] / 2 - 2.0) * st["planes"]["chord"]) if "planes" in st else 2 * 30 * 7, plane_incidence_deg=0.0, plane_x=60.5 * 5 - st["length_mm"] / 2,
        rudder_deg=0.0, screws=[], trim_moment_gmm=0.0,
    )

# ------------------------------------------------------------------ prostředí a nejistoty
FUELS = {
    # evropský prášek do pečiva (soda + pyrofosfát SAPP): celkem ~70 ml/g, za studena jen část, pomalu
    "powder_cold": dict(yield_ml_per_g=70.0, fast_frac=0.30, tau_fast_s=60.0, tau_slow_s=2500.0, dose_g=1.5),
    # rychlý jednočinný prášek (MCP, např. Rumford): ~2/3 do 2 min
    "powder_fast": dict(yield_ml_per_g=70.0, fast_frac=0.70, tau_fast_s=50.0, tau_slow_s=600.0, dose_g=1.5),
    # šumivá tableta / směs soda + kyselina citronová 3:2 hm.: ~150 ml/g, rozpouštění 1–3 min
    "tablet": dict(yield_ml_per_g=150.0, fast_frac=0.90, tau_fast_s=90.0, tau_slow_s=600.0, dose_g=1.0),
}
def env_default(fuel="powder_cold"):
    f = FUELS[fuel]
    return dict(theta_r_deg=50.0, theta_a_deg=90.0, powder_g=f["dose_g"], yield_ml_per_g=f["yield_ml_per_g"], fast_frac=f["fast_frac"],
                tau_fast_s=f["tau_fast_s"], tau_slow_s=f["tau_slow_s"], pla_density=1.24, tank_depth_m=0.30, water_T_C=22.0, fuel=fuel, reseal_rng=0.5)

def capillary_pressure(d_mm, theta_deg):
    return 2 * SIGMA * math.cos(math.radians(theta_deg)) / (d_mm / 2 * 1e-3)

def eotvos(d_mm):
    return RHO * G * (d_mm * 1e-3) ** 2 / SIGMA

# ------------------------------------------------------------------ simulace jednoho běhu
def simulate(des, env, t_end=600.0, dt=1e-3, record=True):
    mm = 1e-3
    A_cup = des["cup_w"] * des["cup_l"] * mm * mm                     # půdorys komory
    V_cap_full = A_cup * (des["cup_top"] - des["plug_top"]) * mm       # plyn, až k zátce
    r_ch = des["siphon_d"] / 2 * mm
    Pcap_ch = capillary_pressure(des["siphon_d"], env["theta_r_deg"])  # plyn vytlačuje vodu z kanálku
    Pcap_hole = capillary_pressure(des["hole_d"], env["theta_r_deg"])
    P_hyst = 2 * SIGMA * (math.cos(math.radians(env["theta_r_deg"])) - math.cos(math.radians(env["theta_a_deg"]))) / r_ch
    Eo = eotvos(des["siphon_d"])
    # hmotnosti a objemy
    m_print = des["mass_print_g"] * 1e-3 * env["pla_density"] / 1.24
    m_screws = sum(s["mass_g"] for s in des["screws"]) * 1e-3
    m_powder = env["powder_g"] * 1e-3
    V_disp = des["disp_cm3"] * 1e-6 + sum(s["mass_g"] for s in des["screws"]) * 1e-3 / 7800 + m_powder / 1600
    m = m_print + m_screws + m_powder
    W = m * G
    # plocha a odpor pro svislý pohyb (trup bokem), přidaná hmotnost válce
    A_plan = des["L_hull"] * des["B_hull"] * mm * mm * 0.85
    Cd_v = 1.1
    m_added_v = RHO * math.pi * (des["B_hull"] / 2 * mm) ** 2 * des["L_hull"] * mm
    A_front = math.pi * (des["B_hull"] / 2 * mm) ** 2 * 0.9
    Cd_x = 0.35; m_added_x = 0.05 * m
    # křídla (hloubková kormidla)
    S_w = des["plane_area_mm2"] * mm * mm; phi = math.radians(des["plane_incidence_deg"])
    # vratný moment sklonu: těžiště pod středem vztlaku
    BG = (des["cb_z"] - des["cg_z"]) * mm
    # kinetika plynu (mol -> objem při tlaku v hloubce)
    V_total_gas_STP = env["powder_g"] * env["yield_ml_per_g"] * 1e-6
    fast_left = V_total_gas_STP * env["fast_frac"]; slow_left = V_total_gas_STP * (1 - env["fast_frac"])
    # stav
    z = -env["tank_depth_m"] + 0.02   # hloubka kýlu (0 = hladina), start u dna
    vz = 0.0; x = 0.0; vx = 0.0; theta = 0.0; psi = 0.0
    V_gas = 0.0; gas_in_channel = False; siphon_primed = True; dumping = False; dump_t = 0.0
    n_cycles = 0; t_last_trig = None; cycle_times = []; events = []
    leak_vent = 0.0; surfaced_time = 0.0; t_first_rise = None; dist = 0.0
    rec = {k: [] for k in ("t", "z", "Vgas", "theta", "x", "vz", "vx", "state")} if record else None
    state = "fill"
    n = int(t_end / dt)
    for i in range(n):
        t = i * dt
        depth_keel = -z
        # geometrie plynu v komoře
        h_gas = V_gas / A_cup                                        # tloušťka vrstvy plynu
        z_if = des["cup_top"] * mm - h_gas                           # hladina v komoře nad kýlem (m)
        wet = z_if > des["plug_top"] * mm + 0.5 * mm                  # prášek je pod vodou?
        # produkce plynu (jen když je prášek namočený); objem při místním tlaku
        P_abs = 101325 + RHO * G * max(0.0, depth_keel - z_if)
        rate = 0.0
        if wet or V_gas < 1e-9:
            r1 = fast_left / env["tau_fast_s"]; r2 = slow_left / env["tau_slow_s"]
            fast_left -= r1 * dt; slow_left -= r2 * dt; rate = (r1 + r2) * 101325 / P_abs
        dV = rate * dt
        # --- sifon: vstup plynu do kanálku a spuštění
        head_gas = RHO * G * h_gas                                   # tlak, kterým plyn tlačí na vodu u stropu
        if not dumping:
            if siphon_primed:
                if head_gas > Pcap_ch: gas_in_channel = True
                # čelo plynu v kanálku je o kapilární výšku výš než hladina v komoře
                z_front = z_if + Pcap_ch / (RHO * G)
                if gas_in_channel and z_front <= des["bend_z"] * mm:
                    dumping = True; dump_t = 0.0; state = "dump"; n_cycles += 1
                    if t_last_trig is not None: cycle_times.append(t - t_last_trig)
                    t_last_trig = t; events.append(("trigger", t, V_gas))
            else:
                # ohyb není zalitý vodou: plyn protéká rovnou ven, komora se neplní
                leak_vent += dV; dV = 0.0; state = "vent"
            # únik dírkami zátky, když je hladina u zátky a plyn tlačí víc než kapilára dírky
            if z_if <= des["plug_top"] * mm + 1e-6 and not dumping:
                over = RHO * G * (des["plug_top"] * mm - z_if) + 0.0
                # přebytek produkce nemá kam: jde dírkami (ztráta) – komora je „plná“
                leak_vent += dV; dV = 0.0; state = "full-leak"
        V_gas += dV
        if dumping:
            # vodní sloupec ve vzestupné větvi se vytlačí, pak plyn uniká; tok omezen odporem kanálku
            dump_t += dt
            dp = RHO * G * (des["exit_z"] * mm - z_if) - Pcap_ch
            if dp > 0:
                A_ch = math.pi * r_ch ** 2
                Q = 0.6 * A_ch * math.sqrt(2 * dp / RHO_GAS)         # výtok plynu clonou
                Q = min(Q, V_gas / dt)
                V_gas -= Q * dt
            if V_gas <= 1e-9 or dump_t > 3.0:
                V_gas = max(V_gas, 0.0); dumping = False; gas_in_channel = False; state = "fill"
                # znovuzalití ohybu: bubliny musí z kanálku ven
                dp_reseal = RHO * G * (des["exit_z"] - des["inlet_z"]) * mm
                p_grav = min(1.0, max(0.0, (Eo - 3.37) / (10.0 - 3.37)))
                siphon_primed = (dp_reseal > P_hyst) or (Eo >= 10.0) or (env.get("reseal_rng", 0.5) < p_grav)
                events.append(("dump_end", t, siphon_primed))
        # --- svislý pohyb
        B = RHO * G * (V_disp + V_gas)
        # u hladiny: vynořená část ubírá vztlak (paluba nad vodou)
        if z > -des["H_deck"] * mm * 0.72:
            emerged = (z + des["H_deck"] * mm * 0.72) * (des["L_hull"] * des["B_hull"] * 0.55 * mm * mm)
            B -= RHO * G * min(emerged, V_disp * 0.3)
        Fz = B - W - 0.5 * RHO * Cd_v * A_plan * vz * abs(vz)
        # křídla: relativní proud (vx, vz) v souřadnicích země, tělo nakloněné o theta (kladné = příď nahoru)
        u_rel = math.hypot(vx, vz)
        Fx_w = Fz_w = 0.0
        if u_rel > 1e-4 and S_w > 0:
            flow_ang = math.atan2(-vz, -vx)             # směr nabíhajícího proudu
            # úhel náběhu desky: osa tětivy = theta + phi
            chord_ang = theta + phi
            alpha = (chord_ang - math.atan2(vz, vx) + math.pi) % (2 * math.pi) - math.pi
            alpha = (alpha + math.pi) % (2 * math.pi) - math.pi
            # deska: CL = 2π sinα·cosα (omezeno), CD = 1.28 sin²α + 0.02
            CL = max(-1.1, min(1.1, 2 * math.pi * math.sin(alpha) * math.cos(alpha)))
            CD = 1.28 * math.sin(alpha) ** 2 + 0.02
            q = 0.5 * RHO * u_rel ** 2 * S_w
            # vztlak kolmo na proud, odpor proti pohybu
            ex, ez = vx / u_rel, vz / u_rel
            Fx_w = q * (-CD * ex + CL * (-ez)); Fz_w = q * (-CD * ez + CL * (ex))
        Fx = Fx_w - 0.5 * RHO * Cd_x * A_front * vx * abs(vx)
        Fz += Fz_w
        az = Fz / (m + m_added_v); ax = Fx / (m + m_added_x)
        vz += az * dt; z += vz * dt; vx += ax * dt; x += vx * dt; dist += abs(vx) * dt
        if z >= 0.0: z = 0.0; vz = min(vz, 0.0); surfaced_time += dt
        floor = -env["tank_depth_m"] + 0.02
        if z <= floor: z = floor; vz = max(vz, 0.0)
        if t_first_rise is None and vz > 0.01 and z > floor + 0.005: t_first_rise = t
        # --- sklon (kvazistaticky): moment plynu + trim vs. vratný moment
        M_gas = RHO * G * V_gas * des["cup_x_offset"] * mm
        M_trim = des["trim_moment_gmm"] * 1e-3 * G * mm
        M_rest = m * G * BG
        theta = max(-0.6, min(0.6, (M_gas + M_trim) / max(M_rest, 1e-9)))
        if record and i % 20 == 0:
            rec["t"].append(t); rec["z"].append(z); rec["Vgas"].append(V_gas * 1e6); rec["theta"].append(math.degrees(theta)); rec["x"].append(x); rec["vz"].append(vz); rec["vx"].append(vx); rec["state"].append(state)
        if fast_left + slow_left < 0.02 * V_total_gas_STP and V_gas < 1e-8 and state == "fill" and t > 30:
            break
    # kormidlo: odhad poloměru zatáčky (kluz) – Nomotův odhad R ≈ 1,6·L / tan(δ)
    R_turn = (1.6 * des["L_hull"] * mm / math.tan(math.radians(max(1e-3, abs(des["rudder_deg"]))))) if des["rudder_deg"] else float("inf")
    fail = None
    if n_cycles == 0: fail = "nespustil se sifon" if leak_vent > 0 else "nestoupla"
    elif not siphon_primed: fail = "sifon se po vypuštění nezalil (trvalé větrání)"
    summary = dict(design=des["name"], cycles=n_cycles, first_rise_s=t_first_rise, mean_cycle_s=(float(np.mean(cycle_times)) if cycle_times else None),
                   surfaced_frac=surfaced_time / max(t, 1e-9), leak_cm3=leak_vent * 1e6, Eo=Eo, Pcap_ch_mmH2O=Pcap_ch / (RHO * G) * 1e3,
                   trigger_margin_mm=des["bend_z"] - Pcap_ch / (RHO * G) * 1e3 - des["plug_top"], net_g=(W - RHO * G * V_disp) / G * 1e3,
                   distance_m=dist, turn_radius_m=R_turn, fail=fail, t_end=t)
    return summary, rec

# ------------------------------------------------------------------ Monte Carlo
def monte_carlo(des, n=300, seed=1, t_end=400.0, fuel="powder_cold"):
    rng = np.random.default_rng(seed); ok = 0; fails = {}; cycles = []; periods = []
    for _ in range(n):
        env = env_default(fuel)
        env["theta_r_deg"] = rng.uniform(35, 60); env["theta_a_deg"] = min(110, env["theta_r_deg"] + rng.uniform(25, 50))
        env["powder_g"] *= rng.uniform(0.7, 1.4); env["yield_ml_per_g"] *= rng.uniform(0.8, 1.2); env["fast_frac"] = min(0.95, env["fast_frac"] * rng.uniform(0.7, 1.3))
        env["pla_density"] = rng.uniform(1.21, 1.27); env["reseal_rng"] = rng.uniform(0, 1)
        d = dict(des); d["screws"] = [dict(mass_g=1.9)] * int(rng.integers(0, 3))
        s, _ = simulate(d, env, t_end=t_end, dt=5e-3, record=False)
        good = s["cycles"] >= 3 and s["fail"] is None
        ok += good
        if not good: fails[s["fail"] or "málo cyklů"] = fails.get(s["fail"] or "málo cyklů", 0) + 1
        else: cycles.append(s["cycles"]); periods.append(s["mean_cycle_s"] or 0)
    return dict(design=des["name"], p_ok=ok / n, fails=fails, cycles_med=(float(np.median(cycles)) if cycles else 0), period_med=(float(np.median(periods)) if periods else 0))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--mc", action="store_true"); ap.add_argument("--plot", action="store_true"); a = ap.parse_args()
    base = design_v2(); env = env_default()
    variants = []
    for d_mm in (3.2, 6.0, 8.0, 9.0, 10.0):
        for bend, top in ((5.5, 12.0), (7.0, 14.0), (8.0, 15.0)):
            v = dict(base); v["name"] = f"d{d_mm}_bend{bend}"; v["siphon_d"] = d_mm; v["bend_z"] = bend; v["cup_top"] = top; v["inlet_z"] = top; v["exit_z"] = top + 3.0 + d_mm / 2 + 1.0; variants.append(v)
    if a.mc:
        for fuel in ("powder_cold", "tablet"):
            print(f"\n=== palivo: {fuel} ===")
            rows = [monte_carlo(v, n=60, fuel=fuel, t_end=240.0) for v in variants]
            print(f"{'návrh':>14} {'P(funguje)':>11} {'cyklů':>6} {'perioda s':>10}  selhání")
            for r in rows: print(f"{r['design']:>14} {r['p_ok']*100:>10.0f}% {r['cycles_med']:>6.0f} {r['period_med']:>10.0f}  {r['fails']}")
            json.dump(rows, open(os.path.join(OUT, f"mc_{fuel}.json"), "w"), indent=1, ensure_ascii=False)
    else:
        s, rec = simulate(base, env, t_end=300.0)
        print(json.dumps(s, indent=1, ensure_ascii=False))
        if a.plot:
            import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
            fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
            ax[0].plot(rec["t"], np.array(rec["z"]) * 100); ax[0].set_ylabel("hloubka kýlu [cm]")
            ax[1].plot(rec["t"], rec["Vgas"]); ax[1].set_ylabel("plyn v komoře [cm³]")
            ax[2].plot(rec["t"], rec["theta"]); ax[2].set_ylabel("sklon [°]"); ax[2].set_xlabel("čas [s]")
            plt.tight_layout(); plt.savefig(os.path.join(OUT, "run_v2.png"), dpi=110)
