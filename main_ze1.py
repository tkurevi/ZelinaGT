"""
main_ze1.py
===========

Nodal-analysis driver for the Zelina-1 (Ze-1) geothermal well and its
planned 600 m doublet partner, Sava depression, Pannonian basin, Croatia.

Ze-1 (drilled 1966, TD 1315 m, Sveti Ivan Zelina spa) produces low-
enthalpy geothermal water (~45-50 C, ~0.87 g/L NaCl, no hydrocarbons)
from the IVA pjescenjaci member of the Ivanic Grad Fm.  The producing
reservoir is SEVEN stacked sands between 753 and 891 m:

    S1 753.0-764.0   S2 773.0-786.0   S3 791.5-800.0   S4 814.0-818.5
    S5 839.0-851.0   S6 853.5-868.5   S7 877.0-891.0      (net 78.0 m)

All seven were logged (1966 suite: sonic, dual normal R16/R64, laterolog,
microlog, SP, caliper, temperature -- NO gamma/neutron/density) and the
reservoir is genuinely MULTILAYER, so unlike OsGT-1 (single sandstone),
Zelina-1 uses the full CommingledReservoir machinery (cf. VGGT-1 2-layer).

-----------------------------------------------------------------------
RESERVOIR PROPERTIES -- all VALIDATED / RE-DERIVED from the raw data
-----------------------------------------------------------------------
* Porosity: sonic (Raymer-Hunt-Gardner, matrix 55.5 us/ft) calibrated
  to core (28-34 % in the 773-789 m cores) -> C=0.601 -> phi 0.297-0.322
  (uniform across all 7 sands).  Re-derived directly from the LAS files.
* Permeability: porosity/Kozeny-Carman distribution anchored to the 1980
  RLT pressure-buildup transmissivity (Option B, 93.9 -> 68.4 mD), i.e.
  the geophysics-distributed FORMATION k that also seeds FEFLOW.
      net kh(Option B) = 6466 mD.m ; T = 1.04e-4 m2/s.
  Independent re-digitisation of the 1980 buildup (Tabela 2) gives a
  middle-time radial-flow slope ~6.9 m/cycle -> T ~9e-5, kh ~5540 mD.m;
  the late slope steepens (partial flow boundary at ~50-90 m, Bourdet/
  Horner).  Honest formation-kh range is 5500-6500 mD.m.
* DELIVERABILITY CALIBRATION (the key well-model choice): the formation
  buildup kh OVER-predicts the OBSERVED 1980 step-rate deliverability
  (~2.5-2.9 m3/h/at, AOF ~777 m3/d at 12.7 at) by ~1.5-2x -- the classic
  transient-kh vs pseudo-steady-deliverability gap.  We therefore keep
  the geophysics-distributed formation k (so the well model and FEFLOW
  share one permeability field) and absorb the gap into a single
  effective deliverability skin, calibrated so the commingled IPR
  reproduces the field PI:

      skin_total = +2.47  ->  J = 2.70 m3/h/at  (matches 1980 step-rate)

  This mirrors the OsGT-1 philosophy (skin calibrated to the field ESP
  test).  s=+2.47 is an EFFECTIVE pseudo-skin standing in for completion/
  partial-penetration/near-well-heterogeneity + the buildup-vs-IPR gap,
  NOT a measured mechanical skin (the DST PR=1.77 implied s~-2.4 locally).
* P_res(z): static gradient (Tabela 3), fit P[at] = 0.1005*z + 6.5,
  which simultaneously reproduces the 6.6 at static WELLHEAD pressure
  (artesian) and the report's ~91 at reservoir pressure at 827 m.
* T_res(z): static gradient (Tabela 3), fit T[C] = 0.0335*z + 23.25
  (48.7 C at S1 top -> 52.9 C at S7).  Produced wellhead ~45-50 C.

-----------------------------------------------------------------------
WELL CONFIGURATION -- decided by the IPR/ESP-depth check (see below)
-----------------------------------------------------------------------
Ze-1 construction: 9-5/8" casing 0-268 m; 7-5/8" hole to TD; cemented
hanging slotted liner 5-1/2" to 1000 m.  An ESP can therefore only sit
inside the 9-5/8" casing, i.e. <= 268 m TVD.

At the design rate 25.9 L/s (~2237 m3/d ~ 3x anything Ze-1 was tested at)
the drawdown is ~34 bar and the dynamic fluid level falls to ~266 m --
right at the 268 m pump seat.  Max sustainable from Ze-1 with the pump at
268 m is only ~24-25 L/s (optimistic linear IPR; worse with non-Darcy).
=> Ze-1 is MARGINALLY infeasible as the 25.9 L/s producer.

Decision (user-confirmed contingency): use Ze-1 for INJECTION and a NEW
production well (9-5/8" casing to 753 m + slotted liner through the
753-891 m reservoir) with the ESP seatable at ~700 m, which delivers
25.9 L/s with ~20 bar intake margin (capacity ~57 L/s).

Configs provided
----------------
* NEWWELL_PRODUCER  -- PRIMARY production config: new well, 9-5/8" to
                       753 m, ESP at 700 m.  Feeds the techno-economics.
* ZE1_PRODUCER_268  -- Ze-1 as producer, ESP at 268 m (the marginal case;
                       kept so the feasibility limit is reproducible).
The reservoir (7 sands, calibrated IPR) is identical in both -- only the
wellbore geometry and ESP depth differ.
"""
from __future__ import annotations
import os, sys, json, copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from pvt import (ppm_to_molality, bar_to_Pa, Pa_to_bar,
                 m3h_to_ls, m3h_to_m3d, fmt_q)
from vlp import WellGeometry, ThermalContext, FluidStream
from ipr_multilayer import Layer, CommingledReservoir
from nodal_multi import (solve_operating_point_multi,
                          print_operating_point_multi,
                          export_summary_csv_multi,
                          export_layer_breakdown_csv,
                          ipr_curve_commingled,
                          vlp_curve_at_depth)

# =====================================================================
# Validated 7-layer IVA reservoir (shared by all configs)
# =====================================================================
NaCl_PPM   = 870.0       # 0.87 g/L NaCl (very fresh)
R_W        = 0.108       # 8-1/2" bit class; gauge hole in the sands
R_E        = 300.0       # drainage radius (default); partial boundary 50-90 m noted
SKIN_DELIV = 2.47        # EFFECTIVE deliverability skin -> J = 2.70 m3/h/at (1980 step-rate)

# top, base, h_net, phi(sonic, core-cal), k_md(OptionB RLT), P_res(bar), T_res(C)
_SANDS = [
    ('S1', 753.0, 764.0, 11.0, 0.322, 93.9, 81.13, 48.7),
    ('S2', 773.0, 786.0, 13.0, 0.319, 90.8, 83.20, 49.4),
    ('S3', 791.5, 800.0,  8.5, 0.317, 88.3, 84.80, 49.9),
    ('S4', 814.0, 818.5,  4.5, 0.313, 84.0, 86.82, 50.6),
    ('S5', 839.0, 851.0, 12.0, 0.309, 79.3, 89.65, 51.6),
    ('S6', 853.5, 868.5, 15.0, 0.310, 81.0, 91.22, 52.1),
    ('S7', 877.0, 891.0, 14.0, 0.297, 68.4, 93.50, 52.9),
]

def reservoir_layers():
    """Build the 7 validated Ze-1 IVA sand Layers (list of dicts for Layer)."""
    return [
        dict(name=nm, top_depth_m=t, bottom_depth_m=b, h_net_m=h,
             k_md=k, P_res_bar=P, T_res_C=T,
             r_w=R_W, r_e=R_E, NaCl_ppm=NaCl_PPM,
             mode='darcy_radial', regime='pss',
             skin_total=SKIN_DELIV, porosity=phi, D_nonDarcy=0.0)
        for (nm, t, b, h, phi, k, P, T) in _SANDS
    ]

# Wellbore static-column density that reproduces the observed 0.1005 at/m
# gradient and the 6.6 at static WHP (cooler upper column than reservoir T).
RHO_WB = 1000.0

# =====================================================================
# NEWWELL_PRODUCER -- PRIMARY config (new production well + ESP @700 m)
# =====================================================================
NEWWELL_PRODUCER = dict(
    case_name='Zelina doublet - NEW production well (9-5/8" to 753 m, ESP @700 m)',
    layers=reservoir_layers(),
    commingled=dict(
        reference_depth_m=None,          # default = deepest sand mid-perf (884 m)
        wellbore_density_kg_m3=RHO_WB,
    ),
    # New well: 9-5/8" casing 0-753 m, slotted liner 753-891 (reservoir).
    # ESP at 700 m inside the 9-5/8" casing; 5-1/2" production tubing above.
    well=dict(
        tubing_ID=0.1240,  tubing_OD=0.1397, roughness=46e-6,   # 5-1/2" tubing
        segments=[
            dict(z_top=0.0,   z_bot=700.0, ID=0.1240, OD=0.1397, roughness=46e-6),  # 5-1/2" tbg above ESP
            dict(z_top=700.0, z_bot=753.0, ID=0.2205, OD=0.2445, roughness=46e-6),  # 9-5/8" casing below ESP
            dict(z_top=753.0, z_bot=884.0, ID=0.1240, OD=0.1397, roughness=46e-6),  # 5-1/2" slotted liner thru reservoir
        ],
    ),
    thermal=dict(
        T_surface_C=12.0,            # Zagreb-region annual avg
        geo_gradient_K_m=0.043,      # Tabela 3 static gradient
        T_BH_C=52.9,                 # formation T at z_ref (884 m)
        U_overall=20.0,
        t_prod_days=1.0,
    ),
    fluid=dict(GWR_std=0.0, NaCl_ppm=NaCl_PPM, w_CO2=0.0),  # ~single-phase fresh water
    operating=dict(
        WHP_bar=3.0,                 # production WHP (gauge+atm); design pumped
        n_segments=50,
        pump=dict(z_intake_m=700.0, dP_bar=32.0),   # ESP duty -> design 25.9 L/s
    ),
)

# =====================================================================
# ZE1_PRODUCER_268 -- Ze-1 as producer (ESP @268 m) : marginal case
# =====================================================================
ZE1_PRODUCER_268 = copy.deepcopy(NEWWELL_PRODUCER)
ZE1_PRODUCER_268['case_name'] = (
    'Zelina - Ze-1 as producer (ESP @268 m, 9-5/8" 0-268 m) : MARGINAL at 25.9 L/s')
# Ze-1 geometry: 9-5/8" 0-268 m, then 7-5/8" hole / 5-1/2" slotted liner to 1000 m.
# ESP at 268 m (deepest possible in 9-5/8"); 4-1/2" tubing above.
ZE1_PRODUCER_268['well'] = dict(
    tubing_ID=0.1005, tubing_OD=0.1143, roughness=46e-6,       # 4-1/2" tubing
    segments=[
        dict(z_top=0.0,   z_bot=268.0, ID=0.1005, OD=0.1143, roughness=46e-6),  # 4-1/2" tbg above ESP
        dict(z_top=268.0, z_bot=753.0, ID=0.1900, OD=0.1937, roughness=46e-6),  # 7-5/8" hole below ESP
        dict(z_top=753.0, z_bot=884.0, ID=0.1240, OD=0.1397, roughness=46e-6),  # 5-1/2" slotted liner thru reservoir
    ],
)
ZE1_PRODUCER_268['operating']['pump'] = dict(z_intake_m=268.0, dP_bar=24.0)

# =====================================================================
# Builders (identical pattern to OsGT-1)
# =====================================================================
def build_commingled(config):
    layers = [Layer(**L) for L in config['layers']]
    return CommingledReservoir(
        layers,
        reference_depth_m=config['commingled'].get('reference_depth_m'),
        wellbore_density_kg_m3=config['commingled'].get('wellbore_density_kg_m3'))

def build_well_to_ref(config, z_ref):
    w = config['well']
    return WellGeometry(
        depth_TVD=float(z_ref),
        tubing_ID=w['tubing_ID'], tubing_OD=w.get('tubing_OD'),
        roughness=w.get('roughness', 46e-6), segments=w.get('segments'))

def build_thermal(config, z_ref):
    th = config['thermal']
    T_surf_K = 273.15 + th['T_surface_C']
    T_BH_K = (273.15 + th['T_BH_C']) if 'T_BH_C' in th \
             else T_surf_K + th['geo_gradient_K_m'] * z_ref
    return ThermalContext(
        T_surface=T_surf_K, geo_gradient=th['geo_gradient_K_m'], T_BH=T_BH_K,
        U_overall=th.get('U_overall', 20.0), k_formation=th.get('k_earth', 2.5),
        alpha_formation=th.get('alpha_earth', 1.0e-6),
        time_seconds=th.get('t_prod_days', 1.0) * 86400.0)

# Default config used by the techno-economic app
CALIBRATED_CONFIG = NEWWELL_PRODUCER


# =====================================================================
# Standalone: print reservoir + nodal operating point (text only)
# =====================================================================
def run_full_analysis(config, out_dir, verbose=True):
    os.makedirs(out_dir, exist_ok=True)
    comm = build_commingled(config)
    well = build_well_to_ref(config, comm.z_ref)
    thermal = build_thermal(config, comm.z_ref)
    GWR = config['fluid']['GWR_std']
    m_NaCl = ppm_to_molality(config['fluid']['NaCl_ppm'])
    WHP = config['operating']['WHP_bar']
    n_seg = config['operating'].get('n_segments', 50)
    pump_cfg = config['operating'].get('pump')
    pump = None
    if pump_cfg is not None:
        pump = dict(z_intake_m=float(pump_cfg['z_intake_m']),
                    dP_Pa=bar_to_Pa(float(pump_cfg['dP_bar'])))
    if verbose:
        print("=" * 70); print(f"CASE: {config['case_name']}"); print("=" * 70)
        print(comm.describe())
        print(f"\nWell: {well}")
        print(f"GWR_std={GWR}, NaCl={config['fluid']['NaCl_ppm']:.0f} ppm, WHP={WHP:.2f} bar")
    with open(os.path.join(out_dir, 'config.json'), 'w') as f:
        json.dump(copy.deepcopy(config), f, indent=2, default=str)
    op = solve_operating_point_multi(
        WHP, well, GWR, m_NaCl, thermal, comm, n_segments=n_seg, pump=pump,
        q_min_m3h=0.5, q_max_m3h=2.0 * comm.AOF() * 3600.0)
    if verbose:
        print_operating_point_multi(op)
    export_summary_csv_multi(op, os.path.join(out_dir, 'summary.csv'))
    export_layer_breakdown_csv(comm, os.path.join(out_dir, 'layer_breakdown.csv'))
    return op, comm, well, thermal


if __name__ == "__main__":
    sel = sys.argv[1] if len(sys.argv) > 1 else 'newwell'
    cfg = ZE1_PRODUCER_268 if sel == 'ze1' else NEWWELL_PRODUCER
    out = os.path.join(HERE, f'ze1_run_{sel}')
    run_full_analysis(cfg, out)
