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
   * @todo Model-specific parameters are added together with the RU/XB dynamics
   * in later increments.
   */
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {}
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
};

#endif
