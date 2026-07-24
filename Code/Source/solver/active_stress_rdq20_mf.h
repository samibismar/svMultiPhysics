// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_RDQ20_MF_H
#define ACTIVE_STRESS_RDQ20_MF_H

#include "active_stress.h"

/**
 * @brief RDQ20-MF mean-field active stress model.
 *
 * @todo Full model documentation (description, equations, references, units)
 * will be added after the complete implementation (increments 3-5).
 */
class RDQ20MF : public ActiveStress {
public:
  /// Model label, used for factory registration and XML selection.
  static inline const std::string label = "RDQ20-MF";

  /// @name State vector layout
  /// @{

  /// Number of regulatory-unit (RU) probability states (entries 0-15).
  static constexpr unsigned int n_ru_states = 16;

  /// Number of crossbridge (XB) moment states (entries 16-19).
  static constexpr unsigned int n_xb_states = 4;

  /// Total number of state variables.
  static constexpr unsigned int n_state_variables = n_ru_states + n_xb_states;

  /**
   * @brief Flat index of the RU probability state P(TL, TC, TR, CC).
   *
   * Each argument is 0 or 1 and denotes the state of, respectively, the left
   * tropomyosin unit, the central tropomyosin unit, the right tropomyosin unit
   * and the central troponin (calcium unbound/bound). The ordering matches the
   * reference implementation's serialization (TL outermost, CC innermost) and
   * spans [0, 15].
   */
  static constexpr int ru_index(int TL, int TC, int TR, int CC) {
    return 8 * TL + 4 * TC + 2 * TR + CC;
  }

  /// Flat index of the XB moment state @p i (in [0, 3]), spanning [16, 19].
  static constexpr int xb_index(int i) {
    return static_cast<int>(n_ru_states) + i;
  }

  /// @}

  /**
   * @brief Model parameters class.
   *
   * Declares the regulatory-unit (RU) and crossbridge (XB) parameters. The
   * active-tension parameters are added in a later increment. The default values
   * are the human body-temperature calibration of the reference implementation,
   * expressed in the reference's own units (see the corresponding members
   * below).
   */
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("Kbasic", 13.0, required);
      add_parameter("Koff", 100.0, required);
      add_parameter("Q", 2.0, required);
      add_parameter("mu", 10.0, required);
      add_parameter("gamma", 12.0, required);
      add_parameter("Kd0", 0.381, required);
      add_parameter("alphaKd", -0.571, required);
      add_parameter("SL0", 2.2, required);

      add_parameter("r0", 134.31, required);
      add_parameter("alpha", 25.184, required);
      add_parameter("mu0_fP", 32.653, required);
      add_parameter("mu1_fP", 0.778, required);
    }
  };

  /**
   * @brief Constructor.
   */
  RDQ20MF() : ActiveStress(n_state_variables) {}

  /**
   * @brief Construct an instance of model parameters.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

protected:
  /**
   * @brief Read model parameters from a parameter object.
   */
  virtual void read_model_specific_parameters(
      const ActiveStressModelParameters &params) override;

  /**
   * @brief Distribute model parameters to all parallel processes.
   */
  virtual void distribute_model_specific_parameters(const CmMod &cm_mod,
                                                    const cmType &cm) override;

  /**
   * @brief Initialize the state vector for a single node.
   *
   * Sets the state to (1, 0, ..., 0), i.e. all probability mass in the RU state
   * P(0, 0, 0, 0) and all crossbridge moments equal to zero.
   *
   * @param[out] state State vector for a single node, to be initialized by
   *   this function.
   */
  virtual void init_local(Vector<double> &state) const override;

  /**
   * @brief Advance in time for a single node.
   *
   * @note Not implemented yet: the RU and XB state updates are added in later
   * increments. Calling this method throws an exception.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override;

  /**
   * @brief Compute the active tension for a single node.
   *
   * @note Temporary: returns zero until the active-tension formula is
   * implemented in a later increment.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state,
                               const double fiber_stretch) const override;

private:
  /// @name Regulatory-unit (RU) dynamics helpers
  /// @{

  /**
   * @brief Compute the central-tropomyosin transition rate for each local RU
   * configuration.
   *
   * Fills @p rates_T, indexed as @c rates_T[TL][TC][TR][CC], where TL, TC and TR
   * are the binary left, central and right tropomyosin states and CC is the
   * central troponin calcium-binding state. Each entry is the rate at which the
   * central tropomyosin changes state for that configuration. Because the rate
   * depends on the neighbour states TL and TR, nearest-neighbour cooperativity
   * is retained through the tracked TL-TC-TR configuration. These rates depend
   * only on the model parameters, not on calcium or stretch.
   */
  void ru_transition_rates_tropomyosin(double (&rates_T)[2][2][2][2]) const;

  /**
   * @brief Advance the 16 RU-state probabilities by one forward-Euler substep.
   *
   * Computes the probability fluxes caused by central-state transitions and the
   * effective boundary-neighbour transitions from the mean-field closure, then
   * updates @p state_RU in place.
   *
   * @param[in] dt Substep size [s].
   * @param[in] rates_T Central-tropomyosin transition rates,
   *   indexed @c rates_T[TL][TC][TR][CC].
   * @param[in] rates_C Troponin transition rates, indexed @c rates_C[CC][TC].
   * @param[in,out] state_RU The 16 RU-state probabilities,
   *   indexed @c state_RU[TL][TC][TR][CC].
   */
  void ru_forward_euler_substep(double dt,
                                const double (&rates_T)[2][2][2][2],
                                const double (&rates_C)[2][2],
                                double (&state_RU)[2][2][2][2]) const;

  /**
   * @brief Advance the four crossbridge moments by one implicit-Euler step.
   *
   * Computes the permissivity and the effective permissive/non-permissive
   * transition rates from the updated RU probabilities, forms the 4x4 linear
   * system for the implicit update and solves it in place for @p state_XB.
   *
   * @param[in] dt Outer time step [s].
   * @param[in] velocity Shortening velocity @f$-\dot{SL}/SL_0@f$ [s^-1].
   * @param[in] rates_T Central-tropomyosin transition rates,
   *   indexed @c rates_T[TL][TC][TR][CC].
   * @param[in] state_RU The updated 16 RU-state probabilities,
   *   indexed @c state_RU[TL][TC][TR][CC].
   * @param[in,out] state_XB The four crossbridge moments, ordered
   *   @f$[\mu_P^0, \mu_P^1, \mu_N^0, \mu_N^1]@f$.
   */
  void xb_implicit_update(double dt, double velocity,
                          const double (&rates_T)[2][2][2][2],
                          const double (&state_RU)[2][2][2][2],
                          double (&state_XB)[4]) const;

  /// @}

  /// @name Interface unit conversions (svMultiPhysics EM units to reference units)
  /// @{

  /// Calcium conversion, millimolar [mM] to micromolar [uM].
  static constexpr double calcium_mM_to_microM = 1.0e3;

  /// Time conversion, milliseconds [ms] to seconds [s].
  static constexpr double time_ms_to_s = 1.0e-3;

  /// @}

  /// RU forward-Euler substep size [s].
  static constexpr double ru_substep = 2.5e-5;

  /// Fixed reference sarcomere length [um] in the length-dependent dissociation
  /// constant (distinct from the parameter SL0).
  static constexpr double kd_reference_sarcomere_length = 2.15;

  /// @name RU model parameters (reference units)
  /// @{

  double Kbasic;  ///< Basic tropomyosin transition rate [s^-1].
  double Koff;    ///< Troponin unbinding rate [s^-1].
  double Q;       ///< Tropomyosin transition-rate asymmetry factor [-].
  double mu;      ///< Calcium-binding cooperativity factor [-].
  double gamma;   ///< Nearest-neighbour cooperativity factor [-].
  double Kd0;     ///< Calcium dissociation constant at reference length [uM].
  double alphaKd; ///< Length dependence of the dissociation constant [uM/um].
  double SL0;     ///< Reference sarcomere length [um]; maps stretch to length.

  double r0;     ///< Combined attachment-detachment rate at zero velocity [s^-1].
  double alpha;  ///< Coefficient of |v| in r(v) = r0 + alpha * |v| [-].
  double mu0_fP; ///< Permissive influx into the zeroth-moment crossbridge state [s^-1].
  double mu1_fP; ///< Permissive influx into the first-moment crossbridge state [s^-1].

  /// @}
};

#endif
