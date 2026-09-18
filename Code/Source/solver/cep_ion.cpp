// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "cep_ion.h"

#include "FE/Common/FEException.h"
#include "all_fun.h"
#include "utils.h"

namespace cep_ion {
namespace {

// Only assembled nodal quantities are communicated. Compact model histories
// remain separate, including at nodes belonging to more than one domain.
template <class StateIndex>
Vector<double> assemble_state(const ComMod &com_mod, const eqType &eq,
                              StateIndex state_index) {
  Vector<double> sum(com_mod.tnNo), count(com_mod.tnNo);
  for (const auto &dmn : eq.dmn) {
    if (dmn.phys != consts::EquationType::phys_CEP)
      continue;
    svmp::check_not_null<svmp::FE::NotInitializedException>(
        dmn.cep.ionic_model, "ionic model was not constructed.");
    const auto &model = *dmn.cep.ionic_model;
    const int row = state_index(model);
    if (row < 0)
      continue;
    const auto &nodes = model.get_node_indices();
    for (size_t column = 0; column < nodes.size(); ++column) {
      sum[nodes[column]] += model.state_value(row, column);
      count[nodes[column]] += 1.0;
    }
  }

  if (com_mod.dmnId.size() != 0) {
    all_fun::commu(com_mod, sum);
    all_fun::commu(com_mod, count);
  }
  for (int node = 0; node < com_mod.tnNo; ++node)
    if (!utils::is_zero(count[node]))
      sum[node] /= count[node];
  return sum;
}

} // namespace

void set_voltage(const eqType &eq, const Vector<double> &voltage) {
  for (const auto &dmn : eq.dmn) {
    if (dmn.phys != consts::EquationType::phys_CEP)
      continue;
    svmp::check_not_null<svmp::FE::NotInitializedException>(
        dmn.cep.ionic_model, "ionic model was not constructed.");
    dmn.cep.ionic_model->set_voltage(voltage);
  }
}

Vector<double> assemble_calcium(const ComMod &com_mod, const eqType &eq) {
  return assemble_state(com_mod, eq, [](const IonicModel &model) {
    return model.get_calcium_index();
  });
}

Vector<double> assemble_output(const ComMod &com_mod, const eqType &eq,
                               const std::string &name) {
  return assemble_state(com_mod, eq, [&name](const IonicModel &model) {
    for (const auto &[field, row] : model.get_output_variables())
      if (field == name)
        return row;
    return -1;
  });
}

void cep_init(Simulation *simulation, SolutionStates &solutions) {
  using namespace consts;
  auto &com_mod = simulation->com_mod;
  auto &Yo = solutions.old.get_velocity();

  for (auto &eq : com_mod.eq) {
    if (eq.phys != EquationType::phys_CEP)
      continue;
    for (int domain = 0; domain < eq.nDmn; ++domain) {
      auto &dmn = eq.dmn[domain];
      if (dmn.phys != EquationType::phys_CEP)
        continue;
      std::vector<int> nodes;
      for (int node = 0; node < com_mod.tnNo; ++node) {
        if (!all_fun::is_domain(com_mod, eq, node, EquationType::phys_CEP))
          continue;
        if (com_mod.dmnId.size() != 0) {
          if (dmn.Id >= 0 && !utils::btest(com_mod.dmnId(node), dmn.Id))
            continue;
        } else if (domain != 0) {
          continue;
        }
        nodes.push_back(node);
      }
      svmp::check_not_null<svmp::FE::NotInitializedException>(
          dmn.cep.ionic_model, "ionic model was not constructed.");
      dmn.cep.ionic_model->initialize_state(com_mod.tnNo, nodes);
    }

    const auto voltage = assemble_state(com_mod, eq,
                                        [](const IonicModel &) { return 0; });
    for (int node = 0; node < com_mod.tnNo; ++node)
      Yo(eq.e, node) = voltage[node];
    // The first reaction starts from the averaged initial voltage, even if
    // VTU initialization or boundary conditions subsequently change Yo.
    set_voltage(eq, voltage);
  }
}

void cep_integ(Simulation *simulation, const int iEq, const int iDof,
               SolutionStates &solutions, const Vector<double> &I4f) {
  auto &Yo = solutions.old.get_velocity();
  static bool IPASS = true;
  auto &com_mod = simulation->com_mod;
  auto &eq = com_mod.eq[iEq];

  // Ignore the first pass, retaining the initialized ionic voltage.
  if (IPASS) {
    IPASS = false;
  } else {
    set_voltage(eq, Yo.row(iDof));
  }

  for (auto &dmn : eq.dmn) {
    if (dmn.phys != consts::EquationType::phys_CEP)
      continue;
    auto &cep = dmn.cep;
    svmp::check_not_null<svmp::FE::NotInitializedException>(
        cep.ionic_model, "ionic model was not constructed.");
    cep.ionic_model->advance_time_step(
        cep.odes, cep.imyo, com_mod.time - com_mod.dt, com_mod.dt, cep.dt,
        cep.Ksac, com_mod.x, I4f,
        [&cep](double time, const Vector<double> &x) {
          return cep.stimulus_value(time, x);
        });
  }

  const auto voltage = assemble_state(com_mod, eq,
                                      [](const IonicModel &) { return 0; });
  simulation->cep_mod.calcium = assemble_calcium(com_mod, eq);
  for (int node = 0; node < com_mod.tnNo; ++node)
    Yo(iDof, node) = voltage[node];
}

std::size_t restart_size(const ComMod &com_mod) {
  std::size_t size = 0;
  for (const auto &eq : com_mod.eq) {
    if (eq.phys != consts::EquationType::phys_CEP)
      continue;
    for (const auto &dmn : eq.dmn) {
      if (dmn.phys != consts::EquationType::phys_CEP)
        continue;
      const auto &model = *dmn.cep.ionic_model;
      size += sizeof(double) * model.get_node_indices().size() *
              (model.nX() + model.nG());
    }
  }
  return size;
}

void write_restart(const ComMod &com_mod, std::ostream &stream) {
  for (const auto &eq : com_mod.eq) {
    if (eq.phys != consts::EquationType::phys_CEP)
      continue;
    for (const auto &dmn : eq.dmn)
      if (dmn.phys == consts::EquationType::phys_CEP)
        dmn.cep.ionic_model->write_state(stream);
  }
}

void read_restart(ComMod &com_mod, std::istream &stream) {
  for (auto &eq : com_mod.eq) {
    if (eq.phys != consts::EquationType::phys_CEP)
      continue;
    for (auto &dmn : eq.dmn)
      if (dmn.phys == consts::EquationType::phys_CEP)
        dmn.cep.ionic_model->read_state(stream);
  }
}

} // namespace cep_ion
