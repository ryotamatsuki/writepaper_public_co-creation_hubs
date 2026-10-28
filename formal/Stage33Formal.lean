import Mathlib

/-!
# Stage 33 proof-critical core

This file formalizes algebraic implications of the frozen two-region
public/private innovation-intermediation model. It does NOT derive the
Cournot equilibrium or certify the global continuous-price equilibrium.
Every economic primitive appearing here is explicit and must be mapped
back to the paper's model independently.
-/

namespace Stage33Formal

noncomputable section

/-- Uniform ability: participation mass of types in [ell,1]. -/
def mass (ell : ℝ) : ℝ := 1 - ell

/-- Uniform ability: first ability moment of participating types. -/
def skill (ell : ℝ) : ℝ := (1 - ell ^ 2) / 2

/-- Exact upper-tail moment decomposition; no equilibrium assumptions. -/
theorem skill_factorization (ell : ℝ) :
    skill ell = mass ell * (1 + ell) / 2 := by
  unfold skill mass
  ring

/-- Extra skill intensity relative to random selection of the same mass. -/
theorem selection_wedge (ell : ℝ) :
    skill ell - mass ell / 2 = ell * (1 - ell) / 2 := by
  unfold skill mass
  ring

theorem selection_wedge_strict {ell : ℝ}
    (hzero : 0 < ell) (hone : ell < 1) :
    mass ell / 2 < skill ell := by
  have hp : 0 < ell * (1 - ell) / 2 := by positivity
  have heq := selection_wedge ell
  linarith

/-- Mass of newcomers, and their unnormalized skill, after free access. -/
theorem new_participant_skill (ell : ℝ) :
    (1 / 2 : ℝ) - skill ell = ell ^ 2 / 2 := by
  unfold skill
  ring

/-- A privately priced hub's interior cutoff in the full-coverage other region. -/
def cutoff (U T : ℝ) : ℝ := (1 / 2 : ℝ) - U / (2 * T)

def interiorFee (U T : ℝ) : ℝ := (U + T) / 2

/-- Higher skill complementarity raises the private fee on the interior branch. -/
theorem private_fee_increases (U T₁ T₂ : ℝ) (hT : T₁ < T₂) :
    interiorFee U T₁ < interiorFee U T₂ := by
  unfold interiorFee
  linarith

/-- At positive U,T, higher complementarity increases the participation cutoff. -/
theorem cutoff_increases (U T₁ T₂ : ℝ)
    (hU : 0 < U) (hT₁ : 0 < T₁) (hT : T₁ < T₂) :
    cutoff U T₁ < cutoff U T₂ := by
  have hT₂ : 0 < T₂ := lt_trans hT₁ hT
  have hp₁ : 0 < 2 * T₁ := by positivity
  have hp₂ : 0 < 2 * T₂ := by positivity
  have hcross : U * (2 * T₁) < U * (2 * T₂) := by
    have hpos : 0 < U * (T₂ - T₁) := mul_pos hU (sub_pos.mpr hT)
    nlinarith
  have hrat : U / (2 * T₂) < U / (2 * T₁) :=
    (div_lt_div_iff₀ hp₂ hp₁).2 hcross
  unfold cutoff
  linarith

theorem participation_decreases {ell₁ ell₂ : ℝ} (h : ell₁ < ell₂) :
    mass ell₂ < mass ell₁ := by
  unfold mass
  linarith

theorem aggregate_skill_decreases {ell₁ ell₂ : ℝ}
    (hzero : 0 ≤ ell₁) (h : ell₁ < ell₂) :
    skill ell₂ < skill ell₁ := by
  have hp : 0 < (ell₂ - ell₁) * (ell₂ + ell₁) / 2 := by
    apply div_pos
    · exact mul_pos (sub_pos.mpr h) (by linarith)
    · norm_num
  have hid : skill ell₁ - skill ell₂ =
      (ell₂ - ell₁) * (ell₂ + ell₁) / 2 := by
    unfold skill
    ring
  linarith

/-- Local public-entry fixed-cost threshold, conditional on other-region coverage. -/
def localThreshold (H C B₀ B₁ ell : ℝ) : ℝ :=
  H + C * (B₀ * ell + B₁ * ell ^ 2 / 4)

/-- Restrained planner threshold for the SAME public entry margin. -/
def plannerThreshold (H C A₀ A₁ B₀ B₁ ell : ℝ) : ℝ :=
  H + C * ((A₀ + B₀) * ell + (A₁ + B₁) * ell ^ 2 / 4)

/-- Foreign welfare externality, not an assumed positive spillover. -/
theorem planner_local_wedge (H C A₀ A₁ B₀ B₁ ell : ℝ) :
    plannerThreshold H C A₀ A₁ B₀ B₁ ell -
      localThreshold H C B₀ B₁ ell =
        C * ell * (A₀ + A₁ * ell / 4) := by
  unfold plannerThreshold localThreshold
  ring

theorem positive_foreign_externality
    (H C A₀ A₁ B₀ B₁ ell : ℝ)
    (hC : 0 < C) (hell : 0 < ell)
    (hForeign : 0 < A₀ + A₁ * ell / 4) :
    localThreshold H C B₀ B₁ ell <
      plannerThreshold H C A₀ A₁ B₀ B₁ ell := by
  have hp : 0 < C * ell * (A₀ + A₁ * ell / 4) := by positivity
  have h := planner_local_wedge H C A₀ A₁ B₀ B₁ ell
  linarith

theorem negative_foreign_externality
    (H C A₀ A₁ B₀ B₁ ell : ℝ)
    (hC : 0 < C) (hell : 0 < ell)
    (hForeign : A₀ + A₁ * ell / 4 < 0) :
    plannerThreshold H C A₀ A₁ B₀ B₁ ell <
      localThreshold H C B₀ B₁ ell := by
  have hp : C * ell * (A₀ + A₁ * ell / 4) < 0 := by
    exact mul_neg_of_pos_of_neg (mul_pos hC hell) hForeign
  have h := planner_local_wedge H C A₀ A₁ B₀ B₁ ell
  linarith

/-- Original Stage 33 conditional municipality threshold. -/
def publicValue (N Delta h g lambda ell : ℝ) : ℝ :=
  -(N * Delta) / 8 + h * ell + (N * Delta / 2 + g * lambda) * ell ^ 2

/-- Cubic numerator after substituting the interior private cutoff. -/
def entryPolynomial (N Delta FG h g lambda U T : ℝ) : ℝ :=
  -(N * Delta + 8 * FG) * T ^ 2 +
    4 * h * T * (T - U) +
    (N * Delta + 2 * g * lambda) * (T - U) ^ 2

/-- Exact polynomialization of the municipality's public-entry surplus. -/
theorem public_entry_polynomial_identity
    (N Delta FG h g lambda U T : ℝ) (hT : T ≠ 0) :
    8 * T ^ 2 *
      (publicValue N Delta h g lambda ((T - U) / (2 * T)) - FG) =
    entryPolynomial N Delta FG h g lambda U T := by
  unfold publicValue entryPolynomial
  field_simp
  ring

/-- Positivity of the cubic numerator is equivalent to beneficial local entry. -/
theorem public_entry_iff_polynomial
    (N Delta FG h g lambda U T : ℝ) (hT : 0 < T) :
    FG ≤ publicValue N Delta h g lambda ((T - U) / (2 * T)) ↔
    0 ≤ entryPolynomial N Delta FG h g lambda U T := by
  have hp : 0 < 8 * T ^ 2 := by positivity
  have heq := public_entry_polynomial_identity N Delta FG h g lambda U T (ne_of_gt hT)
  constructor
  · intro hEntry
    have hnonneg :
        0 ≤ 8 * T ^ 2 *
          (publicValue N Delta h g lambda ((T - U) / (2 * T)) - FG) :=
      mul_nonneg (le_of_lt hp) (sub_nonneg.mpr hEntry)
    rw [heq] at hnonneg
    exact hnonneg
  · intro hPsi
    by_contra hfail
    have hneg : publicValue N Delta h g lambda ((T - U) / (2 * T)) - FG < 0 := by
      have hlt : publicValue N Delta h g lambda ((T - U) / (2 * T)) < FG :=
        lt_of_not_ge hfail
      linarith
    have hprod : 8 * T ^ 2 *
        (publicValue N Delta h g lambda ((T - U) / (2 * T)) - FG) < 0 :=
      mul_neg_of_pos_of_neg hp hneg
    rw [heq] at hprod
    linarith

/-- All four entrants have binary actions, so no unilateral deviation iff every
    flip yields a weakly lower payoff, assuming these are the four flip gains. -/
theorem four_player_candidate_NE_iff
    (dGA dPA dGB dPB : ℝ) :
    (dGA ≤ 0 ∧ dPA ≤ 0 ∧ dGB ≤ 0 ∧ dPB ≤ 0) ↔
      ¬ (0 < dGA ∨ 0 < dPA ∨ 0 < dGB ∨ 0 < dPB) := by
  constructor
  · rintro ⟨hA,hP,hB,hQ⟩ h
    rcases h with h | h | h | h <;> linarith
  · intro h
    constructor
    · by_contra hn
      apply h
      left
      exact lt_of_not_ge hn
    constructor
    · by_contra hn
      apply h
      right; left
      exact lt_of_not_ge hn
    constructor
    · by_contra hn
      apply h
      right; right; left
      exact lt_of_not_ge hn
    · by_contra hn
      apply h
      right; right; right
      exact lt_of_not_ge hn

end

end Stage33Formal
