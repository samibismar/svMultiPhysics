// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SET_BC_H 
#define SET_BC_H 

#include "Simulation.h"
#include "SolutionStates.h"
#include "consts.h"
#include "RobinBoundaryCondition.h"

#include <string>

namespace set_bc {

void calc_der_cpl_bc(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

void cplBC_Integ_X(ComMod& com_mod, const CmMod& cm_mod, const bool RCRflag);

void genBC_Integ_X(ComMod& com_mod, const CmMod& cm_mod, const std::string& genFlag);

void rcr_init(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

void RCR_Integ_X(ComMod& com_mod, const CmMod& cm_mod, int istat);

void set_bc_cmm(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);
void set_bc_cmm_l(ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const SolutionStates& solutions);

void set_bc_cpl(ComMod& com_mod, CmMod& cm_mod, const SolutionStates& solutions);

void set_bc_dir(ComMod& com_mod, SolutionStates& solutions);
void set_bc_dir_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa, Array<double>& lA, Array<double>& lY, int lDof);
void set_bc_dir_w(ComMod& com_mod, const SolutionStates& solutions);
void set_bc_dir_wl(ComMod& com_mod, const bcType& lBc, const mshType& lM, const faceType& lFa, const SolutionStates& solutions);

void set_bc_neu(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

/// @brief Prepare old-geometry basal loads before the predictor for com_mod.cTS.
/// Explicit uses prescribed pressure at evaluation_time or accepted svZeroD pressure.
/// Semi-implicit caches unit-pressure loads; fully implicit defers to set_bc_neu.
void prepare_energy_balance(ComMod& com_mod, const CmMod& cm_mod,
    const SolutionStates& solutions, double evaluation_time);

/// @brief Assemble basal loads with a minus sign in the residual and zero tangent.
/// Semi-implicit scales the cache by the current prescribed/coupled pressure.
/// Fully implicit requires set_bc_neu to refresh it using intermediate geometry.
void assemble_energy_balance(ComMod& com_mod, const bcType& bc);

void set_bc_neu_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa, const SolutionStates& solutions);

void set_bc_rbnl(ComMod& com_mod, const faceType& lFa, const RobinBoundaryCondition& robin_bc,
  const SolutionStates& solutions);

void set_bc_trac_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa, const SolutionStates& solutions);

void set_bc_undef_neu(ComMod& com_mod);

void set_bc_undef_neu_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa);

};

#endif

