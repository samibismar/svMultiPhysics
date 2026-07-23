// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_rdq20_mf.h"

#include "Core/Exception.h"

void RDQ20MF::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  // TODO: read model-specific parameters (added with the RU/XB dynamics in
  // later increments). No parameters to read yet.
}

void RDQ20MF::distribute_model_specific_parameters(const CmMod &cm_mod,
                                                   const cmType &cm) {
  // TODO: distribute model-specific parameters (added with the RU/XB dynamics
  // in later increments). No parameters to distribute yet.
}

void RDQ20MF::init_local(Vector<double> &state) const {
  for (unsigned int i = 0; i < n_states; ++i)
    state[i] = 0.0;

  state[ru_index(0, 0, 0, 0)] = 1.0; // == state[0]
}

void RDQ20MF::advance_time_step_local(const double t, const double dt,
                                      const double calcium,
                                      const double fiber_stretch,
                                      const double fiber_stretch_rate,
                                      Vector<double> &state) const {
  svmp::not_implemented(
      "RDQ20-MF active stress dynamics are not implemented yet: the "
      "regulatory-unit and crossbridge state updates are added in a later "
      "increment.");
}

double RDQ20MF::compute_active_tension_local(const Vector<double> &state,
                                             const double fiber_stretch) const {
  // TEMPORARY: returns zero until the active-tension formula is implemented in
  // a later increment.
  return 0.0;
}

REGISTER_ACTIVE_STRESS_MODEL("RDQ20-MF", RDQ20MF);
