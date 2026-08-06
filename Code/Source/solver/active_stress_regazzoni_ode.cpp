// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_regazzoni_ode.h"

#include <algorithm>
#include <cmath>

void RegazzoniODEActiveStress::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  ActiveStressODE::read_model_specific_parameters(params);

  Kbasic = params.get_scalar("Kbasic");
  Koff   = params.get_scalar("Koff");
  Q      = params.get_scalar("Q");
  mu     = params.get_scalar("mu");
  gamma  = params.get_scalar("gamma");
  Kd0     = params.get_scalar("Kd0");
  alphaKd = params.get_scalar("alphaKd");
  if (alphaKd > 0.0)
    svmp::raise<svmp::ParseException>(
        "RegazzoniODEActiveStress: alphaKd must be <= 0.");
  SL0 = params.get_scalar("SL0");
  kd_reference_sarcomere_length =
      params.get_scalar("kd_reference_sarcomere_length");

  r0     = params.get_scalar("r0");
  alpha  = params.get_scalar("alpha");
  mu0_fP = params.get_scalar("mu0_fP");
  mu1_fP = params.get_scalar("mu1_fP");

  LA   = params.get_scalar("LA");
  LM   = params.get_scalar("LM");
  LB   = params.get_scalar("LB");
  a_XB = params.get_scalar("a_XB");
}

void RegazzoniODEActiveStress::distribute_model_specific_parameters(
    const CmMod &cm_mod, const cmType &cm) {
  ActiveStressODE::distribute_model_specific_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &Kbasic);
  cm.bcast(cm_mod, &Koff);
  cm.bcast(cm_mod, &Q);
  cm.bcast(cm_mod, &mu);
  cm.bcast(cm_mod, &gamma);
  cm.bcast(cm_mod, &Kd0);
  cm.bcast(cm_mod, &alphaKd);
  cm.bcast(cm_mod, &SL0);
  cm.bcast(cm_mod, &kd_reference_sarcomere_length);

  cm.bcast(cm_mod, &r0);
  cm.bcast(cm_mod, &alpha);
  cm.bcast(cm_mod, &mu0_fP);
  cm.bcast(cm_mod, &mu1_fP);

  cm.bcast(cm_mod, &LA);
  cm.bcast(cm_mod, &LM);
  cm.bcast(cm_mod, &LB);
  cm.bcast(cm_mod, &a_XB);
}

void RegazzoniODEActiveStress::init_local(Vector<double> &state) const {
  for (unsigned int i = 0; i < n_state_variables; ++i)
    state[i] = 0.0;
  state[ru_index(0, 0, 0, 0)] = 1.0;
}

Vector<double> RegazzoniODEActiveStress::getf(
    const double /*t*/, const Vector<double> &state,
    const double calcium, const double fiber_stretch,
    const double fiber_stretch_rate) const {

  const double sarcomere_length = SL0 * fiber_stretch;
  const double velocity = -fiber_stretch_rate;

  // ── Central-tropomyosin transition rates (parameter-dependent only) ────────
  const RUArray rates_T = ru_transition_rates_tropomyosin();

  // ── Troponin transition rates (calcium- and stretch-dependent) ─────────────
  const double calcium_on_rate =
      Koff /
      (Kd0 - alphaKd * (kd_reference_sarcomere_length - sarcomere_length)) *
      calcium;
  BinaryPairArray rates_C;
  rates_C[0][0] = calcium_on_rate;   // CC=0->1, TC=0
  rates_C[0][1] = calcium_on_rate;   // CC=0->1, TC=1
  rates_C[1][0] = Koff;              // CC=1->0, TC=0
  rates_C[1][1] = Koff / mu;         // CC=1->0, TC=1

  // ── Deserialize RU state from stage state vector ───────────────────────────
  RUArray s;
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC)
          s[TL][TC][TR][CC] = state[ru_index(TL, TC, TR, CC)];

  // ── Central-unit probability fluxes ───────────────────────────────────────
  RUArray flux_TC, flux_CC;
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC) {
          flux_TC[TL][TC][TR][CC] = s[TL][TC][TR][CC] * rates_T[TL][TC][TR][CC];
          flux_CC[TL][TC][TR][CC] = s[TL][TC][TR][CC] * rates_C[CC][TC];
        }

  // ── Mean-field boundary rates from stage state ─────────────────────────────
  BinaryPairArray rate_left, rate_right;
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC) {
      double fsum = 0.0, psum = 0.0;
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC) {
          fsum += flux_TC[TL][TC][TR][CC];
          psum += s[TL][TC][TR][CC];
        }
      rate_left[TL][TC] = (psum > 1.0e-12) ? fsum / psum : 0.0;
    }
  for (int TR = 0; TR < 2; ++TR)
    for (int TC = 0; TC < 2; ++TC) {
      double fsum = 0.0, psum = 0.0;
      for (int TL = 0; TL < 2; ++TL)
        for (int CC = 0; CC < 2; ++CC) {
          fsum += flux_TC[TL][TC][TR][CC];
          psum += s[TL][TC][TR][CC];
        }
      rate_right[TR][TC] = (psum > 1.0e-12) ? fsum / psum : 0.0;
    }

  // ── Boundary fluxes ────────────────────────────────────────────────────────
  RUArray flux_TL, flux_TR;
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC) {
          flux_TL[TL][TC][TR][CC] = s[TL][TC][TR][CC] * rate_right[TC][TL];
          flux_TR[TL][TC][TR][CC] = s[TL][TC][TR][CC] * rate_left[TC][TR];
        }

  // ── XB scalars from stage RU state ────────────────────────────────────────
  double permissivity = 0.0, flux_PN = 0.0, flux_NP = 0.0;
  for (int TL = 0; TL < 2; ++TL)
    for (int TR = 0; TR < 2; ++TR)
      for (int CC = 0; CC < 2; ++CC) {
        permissivity += s[TL][1][TR][CC];
        flux_PN += s[TL][1][TR][CC] * rates_T[TL][1][TR][CC];
        flux_NP += s[TL][0][TR][CC] * rates_T[TL][0][TR][CC];
      }
  const double k_PN =
      (permissivity >= 1.0e-12) ? flux_PN / permissivity : 0.0;
  const double k_NP =
      ((1.0 - permissivity) >= 1.0e-12) ? flux_NP / (1.0 - permissivity) : 0.0;
  const double r = r0 + alpha * std::abs(velocity);
  const double diag_P = r + k_PN;
  const double diag_N = r + k_NP;

  // ── Assemble full 20-state RHS ─────────────────────────────────────────────
  Vector<double> f(n_state_variables);

  // RU derivatives (16 states)
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC)
          f[ru_index(TL, TC, TR, CC)] =
              -flux_TL[TL][TC][TR][CC] + flux_TL[1-TL][TC][TR][CC]
              -flux_TC[TL][TC][TR][CC] + flux_TC[TL][1-TC][TR][CC]
              -flux_TR[TL][TC][TR][CC] + flux_TR[TL][TC][1-TR][CC]
              -flux_CC[TL][TC][TR][CC] + flux_CC[TL][TC][TR][1-CC];

  // XB derivatives (4 states)
  const double mu_P0 = state[xb_index(0)];
  const double mu_P1 = state[xb_index(1)];
  const double mu_N0 = state[xb_index(2)];
  const double mu_N1 = state[xb_index(3)];

  f[xb_index(0)] = -diag_P * mu_P0 + k_NP * mu_N0 + permissivity * mu0_fP;
  f[xb_index(1)] = -velocity * mu_P0 - diag_P * mu_P1 + k_NP * mu_N1
                   + permissivity * mu1_fP;
  f[xb_index(2)] = k_PN * mu_P0 - diag_N * mu_N0;
  f[xb_index(3)] = k_PN * mu_P1 - velocity * mu_N0 - diag_N * mu_N1;

  return f;
}

double RegazzoniODEActiveStress::compute_active_tension_local(
    const Vector<double> &state, const double fiber_stretch) const {
  const double sarcomere_length = SL0 * fiber_stretch;
  return a_XB * (state[xb_index(1)] + state[xb_index(3)]) *
         fraction_single_overlap(sarcomere_length);
}

RegazzoniODEActiveStress::RUArray
RegazzoniODEActiveStress::ru_transition_rates_tropomyosin() const {
  RUArray rates_T;
  for (int TL = 0; TL < 2; ++TL)
    for (int TR = 0; TR < 2; ++TR) {
      const int n = TL + TR;
      const double closing = Kbasic * std::pow(gamma, 2 - n);
      const double opening = Q * Kbasic * std::pow(gamma, n);
      rates_T[TL][1][TR][0] = closing;
      rates_T[TL][1][TR][1] = closing;
      rates_T[TL][0][TR][0] = opening / mu;
      rates_T[TL][0][TR][1] = opening;
    }
  return rates_T;
}

double RegazzoniODEActiveStress::fraction_single_overlap(
    double sarcomere_length) const {
  const double SL = sarcomere_length;
  const double half = (LM - LB) * 0.5;

  if (SL > LA && SL <= LM)
    return (SL - LA) / half;
  if (SL > LM && SL <= 2.0 * LA - LB)
    return (SL + LM - 2.0 * LA) * 0.5 / half;
  if (SL > 2.0 * LA - LB && SL <= 2.0 * LA + LB)
    return 1.0;
  if (SL > 2.0 * LA + LB && SL <= 2.0 * LA + LM)
    return (LM + 2.0 * LA - SL) * 0.5 / half;
  return 0.0;
}

REGISTER_ACTIVE_STRESS_MODEL("RegazzoniRK4", RegazzoniODEActiveStress);
