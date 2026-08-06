// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_REGAZZONI_ODE_H
#define ACTIVE_STRESS_REGAZZONI_ODE_H

#include "active_stress_ode.h"

#include <array>

/**
 * @brief Experimental RK4 variant of the RDQ20-MF active stress model.
 *
 * This class implements the same RDQ20-MF mean-field sarcomere model as
 * @ref RegazzoniActiveStress, but reformulates all 20 state variables (16 RU
 * probabilities + 4 XB moments) as a single explicit ODE system suitable for
 * use with the @ref ActiveStressODE time-integration framework.
 *
 * @note This class is **experimental**. The RU dynamics are stiff: numerical
 * Jacobian analysis at representative resting and activated states gives
 * dominant eigenvalues in the range −4.7 to −6.9 [1/ms], placing the RK4
 * stability limit at approximately 0.4–0.6 ms. The standard svMultiPhysics
 * timestep of 1 ms is outside the RK4 stability region for this model.
 * @ref RegazzoniActiveStress, which uses forward-Euler substepping for the
 * RU dynamics and implicit Euler for the XB moments, is the production
 * integrator. This class exists to quantify the cost of dropping the custom
 * IMEX structure.
 *
 * **State vector layout:** identical to @ref RegazzoniActiveStress —
 * entries 0–15 are RU probabilities P(TL, TC, TR, CC), entries 16–19 are XB
 * moments [μ_P^0, μ_P^1, μ_N^0, μ_N^1].
 */
class RegazzoniODEActiveStress : public ActiveStressODE {
public:
  /// Factory label.
  static inline const std::string label = "RegazzoniRK4";

  /// @name State vector layout (mirrors RegazzoniActiveStress)
  /// @{
  static constexpr unsigned int n_ru_states = 16;
  static constexpr unsigned int n_xb_states = 4;
  static constexpr unsigned int n_state_variables = n_ru_states + n_xb_states;

  static constexpr unsigned int ru_index(unsigned int TL, unsigned int TC,
                                         unsigned int TR, unsigned int CC) {
    return 8 * TL + 4 * TC + 2 * TR + CC;
  }
  static constexpr unsigned int xb_index(unsigned int i) {
    return n_ru_states + i;
  }
  /// @}

  /**
   * @brief Model parameters.
   *
   * Inherits @c ODE_solver from @ref ActiveStressODE::Parameters. All
   * Regazzoni model parameters are required; @c ru_substep is absent because
   * this class does not substep.
   */
  class Parameters : public ActiveStressODE::Parameters {
  public:
    Parameters() : ActiveStressODE::Parameters(label) {
      constexpr bool required = true;

      add_parameter("Kbasic", 0.013, required);
      add_parameter("Koff", 0.1, required);
      add_parameter("Q", 2.0, required);
      add_parameter("mu", 10.0, required);
      add_parameter("gamma", 12.0, required);
      add_parameter("Kd0", 3.81e-4, required);
      add_parameter("alphaKd", -5.71e-4, required);
      add_parameter("SL0", 2.2, required);
      add_parameter("kd_reference_sarcomere_length", 2.15, required);

      add_parameter("r0", 0.13431, required);
      add_parameter("alpha", 25.184, required);
      add_parameter("mu0_fP", 0.032653, required);
      add_parameter("mu1_fP", 7.78e-4, required);

      add_parameter("LA", 1.25, required);
      add_parameter("LM", 1.65, required);
      add_parameter("LB", 0.18, required);
      add_parameter("a_XB", 22.894, required);
    }
  };

  RegazzoniODEActiveStress() : ActiveStressODE(n_state_variables) {}

  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

protected:
  virtual void read_model_specific_parameters(
      const ActiveStressModelParameters &params) override;

  virtual void distribute_model_specific_parameters(const CmMod &cm_mod,
                                                    const cmType &cm) override;

  virtual void init_local(Vector<double> &state) const override;

  /**
   * @brief Compute the 20-state ODE RHS.
   *
   * Returns the derivatives of all 16 RU probabilities and all 4 XB moments
   * evaluated at the current stage state. The mean-field boundary rates
   * (rate_left, rate_right) are recomputed from the stage state at every call,
   * so this function is correct as a stage evaluator for RK4.
   */
  virtual Vector<double> getf(const double t, const Vector<double> &state,
                              const double calcium, const double fiber_stretch,
                              const double fiber_stretch_rate) const override;

  virtual double
  compute_active_tension_local(const Vector<double> &state,
                               const double fiber_stretch) const override;

private:
  using RUArray =
      std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2>;
  using BinaryPairArray = std::array<std::array<double, 2>, 2>;
  using XBArray = std::array<double, 4>;

  RUArray ru_transition_rates_tropomyosin() const;
  double fraction_single_overlap(double sarcomere_length) const;

  /// @name RU parameters
  double Kbasic, Koff, Q, mu, gamma, Kd0, alphaKd, SL0;
  double kd_reference_sarcomere_length;
  /// @}

  /// @name XB parameters
  double r0, alpha, mu0_fP, mu1_fP;
  double LA, LM, LB, a_XB;
  /// @}
};

#endif
