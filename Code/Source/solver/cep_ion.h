// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_ION_H
#define CEP_ION_H

#include "Array.h"
#include "ComMod.h"
#include "Simulation.h"
#include "SolutionStates.h"

#include "all_fun.h"
#include "consts.h"

#include <string>

namespace cep_ion {

/// Initialize each domain's independent state in ascending rank-local node
/// order. Average only voltage, supplying it to the models and old solution.
void cep_init(Simulation *simulation, SolutionStates &solutions);

void cep_integ(Simulation *simulation, const int iEq, const int iDof,
               SolutionStates &solutions, const Vector<double> &I4f);

/// Supply one rank-local voltage field to all participating domain models.
void set_voltage(const eqType &eq, const Vector<double> &voltage);

/// Average each domain's calcium/proxy using its own calcium index.
Vector<double> assemble_calcium(const ComMod &com_mod, const eqType &eq);

/// Average a named ionic output over domains that register it. Nodes without
/// contributors return zero. All ranks must call this, including empty ones.
Vector<double> assemble_output(const ComMod &com_mod, const eqType &eq,
                               const std::string &name);

/// Native binary restart payload: equations and domains in input order, each
/// model's ordinary states followed by gates (see IonicModel::write_state).
/// Sizes come from the initialized layout; the mesh, partition, domain order,
/// and model configuration must match the writing run.
/// No per-model sizes or node IDs are written.
std::size_t restart_size(const ComMod &com_mod);
void write_restart(const ComMod &com_mod, std::ostream &stream);
void read_restart(ComMod &com_mod, std::istream &stream);
}; // namespace cep_ion

#endif
