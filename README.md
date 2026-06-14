# Zelina-1 (Ze-1) geothermal doublet — techno-economic model

A full reservoir → production → surface → economics model for the **Ze-1**
low-enthalpy geothermal well (IVA pješčenjaci, Sava depression, Pannonian
basin), built on the same engine as the OsGT-1 (Osijek) web-app but driven
by **validated Ze-1 data** and a **600 m doublet** with a **new production
well**. Communicate in English; Croatian technical terms are kept.

---

## 1. The well and the reservoir

Ze-1 (drilled 1966, TD 1315 m, Sveti Ivan Zelina spa) produces ~45–53 °C
water (~0.87 g/L NaCl, no hydrocarbons — UV-negative) from **seven stacked
sands, 753–891 m, net 78.0 m**:

| sand | top–base (m) | h_net (m) | φ (–) | k (mD) | P_res (bar) | T_res (°C) |
|------|--------------|-----------|-------|--------|-------------|------------|
| S1 | 753.0–764.0 | 11.0 | 0.322 | 93.9 | 81.1 | 48.7 |
| S2 | 773.0–786.0 | 13.0 | 0.319 | 90.8 | 83.2 | 49.4 |
| S3 | 791.5–800.0 |  8.5 | 0.317 | 88.3 | 84.8 | 49.9 |
| S4 | 814.0–818.5 |  4.5 | 0.313 | 84.0 | 86.8 | 50.6 |
| S5 | 839.0–851.0 | 12.0 | 0.309 | 79.3 | 89.7 | 51.6 |
| S6 | 853.5–868.5 | 15.0 | 0.310 | 81.0 | 91.2 | 52.1 |
| S7 | 877.0–891.0 | 14.0 | 0.297 | 68.4 | 93.5 | 52.9 |

**net kh = 6466 mD·m**, net-average k = 82.9 mD.

All properties were **re-derived from the raw 1966 logs and the 1980 RLT
well test** (not taken on trust):

* **Porosity** — sonic Raymer-Hunt-Gardner (matrix 55.5 µs/ft) calibrated
  to core (28–34 % in the 773–789 m cores) → C = 0.601 → φ 0.297–0.322,
  uniform across the seven sands. (1966 suite has **no** gamma/neutron/
  density, so sonic is the only porosity log.)
* **Permeability** — porosity/Kozeny-Carman distribution **anchored to the
  1980 buildup transmissivity** (Option B). This is the geophysics-
  distributed *formation* k that also seeds FEFLOW. Independent re-
  digitisation of the buildup gives a middle-time radial-flow slope
  ~6.9 m/cycle → kh ~5540 mD·m, with the late slope steepening (partial
  flow boundary at ~50–90 m). Honest formation-kh range **5500–6500 mD·m**.

## 2. The key calibration — deliverability skin

The *formation* buildup kh **over-predicts the observed 1980 step-rate
deliverability** (~2.5–2.9 m³/h/at; AOF ~777 m³/d at 12.7 at) by ~1.5–2×
— the classic transient-kh vs pseudo-steady-deliverability gap. We keep
the geophysics k (so the well model and FEFLOW share one permeability
field) and absorb the gap into a single **effective deliverability skin**:

> **skin = +2.47 → J = 2.70 m³/h/at**, reproducing the field step-rate.

This mirrors the OsGT-1 philosophy (skin calibrated to the field test).
The +2.47 is an *effective* pseudo-skin (completion / partial penetration
/ near-well heterogeneity + the buildup-vs-IPR gap), **not** a measured
mechanical skin (the DST PR = 1.77 implies s ≈ −2.4 locally).

P_res(z) and T_res(z) use the 1980 static gradient (Tabela 3):
`P[at] = 0.1005·z + 6.5` (which also reproduces the **6.6 at artesian
wellhead pressure**) and `T[°C] = 0.0335·z + 23.25`.

## 3. The decisive result — can the ESP sit at 268 m in Ze-1?

Ze-1 construction: **9⅝″ casing to 268 m**, 7⅝″ hole to TD, cemented
slotted 5½″ liner to 1000 m. An ESP can therefore only sit **≤ 268 m**.

At the **design rate 25.9 L/s (≈ 2237 m³/d — about 3× anything Ze-1 was
ever tested at)**:

* drawdown ≈ **34 bar**, flowing BHP ≈ 60 bar @ z_ref (884 m);
* the **natural dynamic fluid level falls to ≈ 286 m** — *below* the 268 m
  pump seat;
* the most Ze-1 can deliver with the pump held at 268 m (≈ 2 bar intake)
  is **≈ 24–25 L/s** on the *optimistic linear* IPR, and **less** once
  non-Darcy/boundary losses at 3× the tested rate are included.

→ **Ze-1 is marginally-to-clearly infeasible as the 25.9 L/s producer.**

**Decision (your confirmed contingency):** use **Ze-1 for injection** and a
**new production well** (9⅝″ casing to 753 m + slotted liner through the
753–891 m reservoir, **drilled to 900 m**). The flowing dynamic level at
25.9 L/s is ~286 m, so the **ESP is set at ~386 m** (≈100 m submergence for
NPSH — *not* all the way down; the static lift the pump does is the dynamic
level itself, so a deeper setting only wastes cable/tubing). Intake ~11 bar,
NPSHa ~110 m. Open-flow capacity ~57 L/s. The wells can be switched.

## 4. The surface plant (single operating point)

HP-assisted **ambient loop**, two operating points (the old Osijek
60/40/15 °C A/B/C are removed):
* **Scenario A (primary, default): 15 °C supply / 10 °C return** → injection 13 °C.
* **Scenario B (optional): 30 °C supply / 20 °C return** → injection 23 °C.

Consumer-side heat pumps **SPF 4.7**; consumer 150 m away.

* **Brine inlet** = wellhead T. The seven sands (48.7–52.9 °C) commingle;
  the bottomhole flowing temperature is the **flow-weighted mixing T, not
  the deepest layer**. At 25.9 L/s the per-layer split (S1→S7) is roughly
  15.4 / 17.8 / 11.4 / 5.8 / 14.9 / 19.2 / 15.3 % → **mixing T = 50.8 °C**,
  cooling only ~0.5 °C up the wellbore (high rate) → **wellhead ≈ 50.3 °C**.
  (The 1980 ~45 °C was at ~7× lower rate, hence more cooling.)
* **Injection (return) temperature is calculated** from the brine being
  cooled to the ambient-loop return plus a cold approach:
  **T_inj ≈ 10 + 3 = 13 °C**. Brine cools 50.3 → 13 °C → **Q_geo ≈ 3.99 MW**
  (HP evaporator duty).
* Consumer heat pumps upgrade it: **Q_delivered = Q_geo·SPF/(SPF−1)
  ≈ 5.07 MW**; **HP electricity = Q_geo/(SPF−1) ≈ 1.08 MW** (with the
  50.3 °C wellhead, Q_geo ≈ 3.99 MW).

The OsGT-1 heat-plant module (which heats a network *toward* brine
temperature, the 2G/3G A/B/C convention) is reused only for PHE/pipeline
hydraulics; the geothermal duty and injection T are set explicitly for the
ambient-loop topology in `economy_Ze1.run_engineering`.

## 5. Doublet

Producer (new) ↔ injector (Ze-1) **spacing 600 m**, circulation 25.9 L/s,
reinjection 13 °C. Layered Gringarten-Sauty decline over the seven sands.

**On the flat decline curve:** the earliest *thermal* breakthrough (most
permeable sand) is at **~22.6 yr** (fluid front ~10 yr × thermal retardation
R_th ≈ 2.2); the cold-front-advance plot now uses that retardation, so the
front reaches the producer at ~yr 23 — consistent with the per-layer
breakthrough bars. **Yet the produced temperature stays ~50.8 °C for >30 yr**
because at 600 m in these thin (4.5–15 m) sands the cap/base rock reheats the
cold front (the Gringarten-Sauty heat-loss term, λ ≈ 0.26). This is real, not
a bug: at 300 m spacing the model drops 8.9 °C in 60 yr, at 150 m it drops
17 °C in 30 yr; an adiabatic 600 m case drops 7 °C. So **600 m is thermally
very safe** — the decline is genuine but slow.

## 6. Files

| file | role |
|------|------|
| `main_ze1.py` | reservoir (7-layer commingled IPR) + wellbore geometries + nodal. Configs `NEWWELL_PRODUCER` (primary) and `ZE1_PRODUCER_268` (the marginal Ze-1-as-producer case). |
| `economy_Ze1.py` | master techno-economic model (`run(Config())`): reservoir bridge, ESP + injection pump, surface plant (HP-assisted), doublet decline, CAPEX/OPEX, NPV/IRR/LCOH. |
| `doublet_decline_Ze1.py` | layered Gringarten-Sauty (+ optional Barends) thermal decline. |
| `geothermal_HE_Ze1.py` | PHE + buried-pipeline hydraulics/heat loss. |
| `esp_geothermal_Ze1.py` | ESP / injection-pump electrical sizing. |
| `vlp.py`, `ipr.py`, `ipr_multilayer.py`, `nodal.py`, `nodal_multi.py`, `pvt.py`, `wellbore_T.py` | reusable well-physics (verbatim from the OsGT engine). |
| `streamlit_app.py` | web UI (Reservoir / Pump / DHS / Economy tabs). Run: `streamlit run streamlit_app.py`. |
| `plotting.py`, `plot_style.py`, `app_figures.py`, `doublet_viz.py`, `geothermo_props.py` | plotting / helpers. |

Run the engine headless: `python main_ze1.py newwell` (or `ze1`), or
`python -c "import economy_Ze1 as E; E.run(E.Config())"`.

## 7. Open assumptions to confirm (the levers that matter)

1. **Heat price / HP ownership — the dominant economic lever.**
   `heat_price_eur_MWhth = 35.0` (your value). Two independent tick-box
   toggles (like `co2_as_revenue`) select the business model and keep it
   internally consistent — the **revenue basis follows who operates the HPs**:
   * `hp_capex_in_scope` (Y/N) — include consumer/central HP CAPEX.
   * `hp_elec_in_opex` (Y/N) — include HP electricity in OPEX **and** sell the
     UPGRADED useful heat; OFF ⇒ plant sells only the LOW-T geothermal heat
     and the consumer owns/operates the HPs.

   At 35 €/MWh the project is uneconomic in every configuration (the cheapest,
   consumer-owned, has **LCOH ≈ 100 €/MWh**), which is the honest signal that a
   single small low-T doublet with a *new* production well needs either a much
   higher heat price, more full-load hours, or a larger/closer heat load:

   | configuration (capex/elec) | delivered MWh/yr | CAPEX | OPEX/yr | LCOH | NPV |
   |---|---|---|---|---|---|
   | Plant owns+operates (ON/ON) | 10 142 | 8.59 M | 541 k | 135 | −12.3 M |
   | Consumer owns+operates (OFF/OFF) | 7 984 | 6.29 M | 245 k | 100 | −5.1 M |
   | Plant owns, consumer pays elec (ON/OFF) | 7 984 | 8.59 M | 250 k | 120 | −7.5 M |
2. **Operating hours** `FLH = 2000 h/yr` and **operating_months = 12**
   (single operating point ⇒ peak = avg = 25.9 L/s). Adjust for the real
   seasonal heat demand.
3. **Ze-1 as injector at 25.9 L/s** (well depth **900 m**). Injection
   BHP = static 93.5 + 34 bar overpressure = 127 bar (injectivity mirrors
   productivity). The injection pump is credited the **gravity assist of the
   cold column (~87 bar)**, so it supplies only **~40 bar at the wellhead**
   (≈177 kW). The artesian overpressure is carried as *resistance* in the
   127 bar BHP — **not** a spurious wellhead bonus. Same ~3×-tested-rate
   challenge in reverse; monitor (the 50–90 m boundary could raise injection
   pressure over time).
7. **ESP setting depth.** Set from the flowing dynamic level (~286 m) + a
   100 m submergence margin → ~386 m (cable/tubing priced to that), not the
   full 700 m. Adjust `ESP_SUBMERGENCE_M` in `get_v3_results` if you want a
   different NPSH margin.
4. **Wellbore cooling.** Ramey conduction on the **50.8 °C mixing T** gives
   ~50.3 °C wellhead at 25.9 L/s (column treated liquid-full; ~0.5 °C loss).
   If the real completion loses more heat, lower the wellhead T (and Q_geo).
5. **New-well cost** `well_cost_eur_per_m = 3950`, depth **900 m** (9 m below
   the S7 base at 891 m) → ~3.56 M€ — the single largest CAPEX item, since the
   new well is forced by the ESP-depth constraint.
6. Rock thermal properties use Pannonian/Osijek defaults; r_w = 0.108 m;
   drainage r_e = 300 m (partial boundary 50–90 m noted separately).
