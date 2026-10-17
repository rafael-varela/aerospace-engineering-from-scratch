# P04.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` (and the probe scripts used to write this file) and given to the figures shown. Every input labelled illustrative is illustrative; nothing here is a measurement.

## 1. Photons and doses ($s = 15$)

(a) n_u = 4 + (15 mod 4) = 7. `transition_wavelength(7, 2)` = **397.12 nm** (vacuum), photon energy **3.122 eV**. Check: E₇ − E₂ = 13.5983(1/4 − 1/49) = 3.1220 eV ✓.

(b) 20 mJ at 532 nm: photon energy 2.3305 eV, so **5.356 × 10¹⁶ photons** (the same count as 10 mJ at 1064 nm, because halving the wavelength doubles the energy per photon). A 15 µs round trip is **2248 m**.

(c) He⁺: −Z² × 13.6057 = **−54.42 eV**.

(d) ΦS = 10¹⁰ × 34.6 = 3.46 × 10¹¹ MeV/g × 1.602 × 10⁻¹⁰ Gy/(MeV/g) = **55.4 Gy = 5.54 krad(Si)**.

**Common mistakes.** Air instead of vacuum wavelengths when comparing with a table; forgetting the factor ½ in the lidar range; mixing g and kg in the dose (1 MeV/g is 1.602 × 10⁻¹⁰ Gy, not 10⁻¹³).

## 2. The chain, its peak and its equilibrium

(a) N₂′ + λ₂N₂ = λ₁N₁(0)e^{−λ₁t}. Multiply by e^{λ₂t}: (N₂e^{λ₂t})′ = λ₁N₁(0)e^{(λ₂−λ₁)t}. Integrate from 0 with N₂(0) = 0: N₂e^{λ₂t} = λ₁N₁(0)(e^{(λ₂−λ₁)t} − 1)/(λ₂ − λ₁). Multiply by e^{−λ₂t}: **N₂ = N₁(0)λ₁/(λ₂ − λ₁)(e^{−λ₁t} − e^{−λ₂t})**. ∎

(b) λ₂N₂ ∝ e^{−λ₁t} − e^{−λ₂t}; its derivative −λ₁e^{−λ₁t} + λ₂e^{−λ₂t} = 0 gives e^{(λ₂−λ₁)t} = λ₂/λ₁, **t_max = ln(λ₂/λ₁)/(λ₂ − λ₁)**. ∎ (For the lesson's half-lives 2 and 0.3, t_max = 0.966.)

(c) λ₂N₂/(λ₁N₁) = λ₂/(λ₂ − λ₁) × (1 − e^{−(λ₂−λ₁)t}) → **λ₂/(λ₂ − λ₁)** as t → ∞ when λ₂ > λ₁; if λ₁ ≪ λ₂ this → 1 (secular equilibrium). Sr-90 → Y-90 (in hours: λ₁ = 2.75 × 10⁻⁶, λ₂ = 1.083 × 10⁻² h⁻¹): **t_max = 765 h ≈ 32 days**, limiting ratio **1.000 25**. A Sr-90 source a few months old has Y-90 activity equal to its own.

(d) One gram holds N_A/M atoms and decays at λN_A/M per second, each releasing Q: **P/m = λN_AQ/M**. Pu-238: λ = 2.505 × 10⁻¹⁰ s⁻¹, N_A/M = 2.530 × 10²¹ g⁻¹, Q = 5.593 MeV = 8.962 × 10⁻¹³ J: **0.568 W/g**.

**Common mistakes.** Losing the sign in λ₂ − λ₁ (the formula holds for either ordering; only the equilibrium needs λ₁ < λ₂); mixing years and hours in (c).

## 3. Americium instead of plutonium ($s = 15$: 25-year mission)

(a) Q = (241.056 829 − 237.048 173 − 4.002 603) u × 931.494 = **5.638 MeV**; specific power **0.1146 W/g**, 4.96 times less than Pu-238 (longer half-life, almost the same Q and M).

(b) Pu-238: 25 years leaves 82.07 % of the heat; inventory **3.93 kg**. Am-241: 96.07 % remains; inventory **16.65 kg**.

(c) Equal inventories require 2^{t/87.7}/2^{t/432.6} = w_Pu/w_Am = 4.955, so t = ln 4.955/(ln 2 (1/87.7 − 1/432.6)) = **254 years** (24.0 kg each). No mission is that long, so for any real mission Pu-238 is lighter by a factor of 4 to 5; Am-241 is attractive for availability, not mass.

**Check.** At 14 years the ratio is 16.36/3.61 = 4.5, smaller than 4.96 because Pu-238 has decayed more; it shrinks to 1 at 254 years.

**Common mistakes.** Using the Pu-238 half-life for both; comparing thermal powers instead of inventories.

## 4. A shielding trade ($s = 15$: 5 krad/yr behind 3 mm; MODEL, illustrative inputs)

(a) L = 0.54/ln(20/5) = **0.3895 g/cm²**, D₀ = **40.0 krad/yr**; half-value thickness L ln 2 = 0.270 g/cm² = **1.00 mm**.

(b) Allowed rate = rating/20.

| rating | areal density | thickness | mass of 600 cm² shell |
|---|---|---|---|
| 30 krad | 1.279 g/cm² | **4.74 mm** | **0.77 kg** |
| 100 krad | 0.810 g/cm² | **3.00 mm** | **0.49 kg** |

(The 100 krad case lands exactly on the second fitted point: 100/20 = 5 krad/yr.)

(c) 3 mm point at 4.0 krad/yr (−20 %): **4.22 mm**; at 6.0 (+20 %): **5.30 mm**. The thickness moves by −11 %/+12 %, *less* than the input, because it depends on the dose through a logarithm and the target is only 1.7 mm beyond the fitted range; a target further out would amplify the uncertainty in L.

(d) The single exponential ignores (i) the hard proton component, which attenuates far more slowly and dominates at depth, and (ii) bremsstrahlung photons generated as electrons stop in the aluminium, which are very penetrating; both make the real dose-depth curve flatten, so the model underestimates the dose behind thick shields.

**Common mistakes.** Using thickness in mm in the exponent with L in g/cm²; forgetting the margin; quoting the model's D₀ as the unshielded dose.

## 5. Upsets and scrubbing ($s = 15$: σ = 2.5 × 10⁻¹⁴ cm²/bit; illustrative)

(a) r = 2.5 × 10⁻¹⁴ × 1 × 4.295 × 10⁹ × 86 400 = **9.28 upsets/day**. One orbit: mean 0.612, P(none) = e^{−0.612} = **0.542**. 100 000 Poisson draws with `default_rng(15)`: fraction with zero events 0.5416, within 0.0006 of the formula (one standard error is 0.0016).

(b) r_w = σφ × 64 = 1.6 × 10⁻¹² s⁻¹ per word; N_w = 6.71 × 10⁷ words.

| scrub interval T | P(some word has two upsets) |
|---|---|
| 1 s | 8.6 × 10⁻¹⁷ |
| 1 min | 3.1 × 10⁻¹³ |
| 1 h | 1.1 × 10⁻⁹ |
| 1 day | 6.4 × 10⁻⁷ |

(c) N_w(r_wT)²/2 = 10⁻⁶ gives **T = 1.08 × 10⁵ s ≈ 30 h**: a daily scrub suffices with margin. The approximation uses P(≥2) ≈ (r_wT)²/2 for one word, accurate when r_wT ≪ 1 (here 1.7 × 10⁻⁷), and adds the words' probabilities, accurate when the total is ≪ 1.

**Common mistakes.** Using the whole-memory rate in the per-word formula; assuming upsets in different words matter (the code corrects each word separately).

## 6. Matching the mitigation to the effect

(a) TID; with margin 2 the design dose is 8 krad against a 10 krad rating: accept it (or add 1 mm of shielding for comfort), the margin is met.
(b) Displacement damage (protons displace silicon atoms in the cells); design margin: size the array for end-of-life power.
(c) SEU in memory; error-correcting code and scrubbing.
(d) SEL; current limiting with automatic power cycling (and a part screened for latch-up).
(e) TID, but the shield was sized with a model used far outside its fitted range: at 12 mm the two-point exponential predicts 0.003 krad/yr, which is meaningless because the real curve flattens; rely on the 300 krad part's hardness and a dose-depth curve for the Jovian environment, not on the model.
(f) Surface or internal charging, outside this lesson's model; follow the charging handbook (`@nasa-hdbk-4002b-in-space-charging`): grounding, conductive coatings, shielding of cables.

## Check-yourself answers

1. E_n = −13.6 eV/n²; H-α (3 → 2) is 656.47 nm in vacuum. The reduced mass lowers the energies by 1 part in 1836 (−13.5983 instead of −13.6057 eV), which moves lines by 0.05 % — 0.36 nm at H-α, far above spectrometer resolution.
2. 5.356 × 10¹⁶ photons (1.165 eV each); 1000 m.
3. B = (Zm_H + Nm_n − M)c²; 8.790 MeV per nucleon for Fe-56. Pu-238 → U-234 + He-4 has Q = +5.593 MeV: the products are more tightly bound (He-4 is exceptionally tight), so the decay releases energy.
4. N = N₀e^{−λt}, T½ = ln 2/λ; after 14 years 2^{−14/87.7} = 89.5 % remains.
5. λN_AQ/M = 2.505 × 10⁻¹⁰ s⁻¹ × 2.530 × 10²¹ g⁻¹ × 8.962 × 10⁻¹³ J = 0.568 W/g.
6. Gray = J/kg of absorbed energy; sievert = gray weighted by the radiation's biological effect (alpha × 20). 1 krad = 10 Gy.
7. Total ionizing dose, displacement damage, single-event effects; shielding helps least against single-event effects from heavy ions.
8. L ln 2 (0.86 mm of aluminium for the lesson's fit). The model's D₀ is an extrapolation of the fitted range back to zero; the soft particles stopped in the first fraction of a millimetre are not in it, so the real unshielded dose is far higher.
