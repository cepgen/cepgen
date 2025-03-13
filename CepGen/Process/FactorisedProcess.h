/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Process_FactorisedProcess_h
#define CepGen_Process_FactorisedProcess_h

#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Process/Process.h"

namespace cepgen::proc {
  /// Generic parton emission-factorised process
  /// \note 0 to 2 dimensions may be used for the scattered diffractive system('s) invariant mass definition.
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2023
  class FactorisedProcess : public Process {
  public:
    /// Class constructor
    /// \param[in] params Parameters list
    /// \param[in] output Produced final state particles
    explicit FactorisedProcess(const ParametersList& params, const spdgids_t& output);
    FactorisedProcess(const FactorisedProcess&);

    double computeWeight() override;
    void fillKinematics() final;

    static ParametersDescription description();

    bool validatedBeamKinematics();

  protected:
    void addEventContent() override;
    void prepareKinematics() final;

    virtual void prepareFactorisedPhaseSpace() = 0;  ///< Prepare central part of the Jacobian after kinematics is set
    virtual double computeFactorisedMatrixElement() = 0;  ///< Factorised matrix element (event weight)

    //--- Mandelstam variables
    double that() const;  ///< \f$\hat t=\frac{1}{2}\left[(p_1-p_3)^2+(p_2-p_4)^2\right]\f$
    double uhat() const;  ///< \f$\hat u=\frac{1}{2}\left[(p_1-p_4)^2+(p_2-p_3)^2\right]\f$

    /// Kinematic variables generator for the phase space coverage
    const std::unique_ptr<PhaseSpaceGenerator> phase_space_generator_;
    const bool symmetrise_;
    const bool store_alphas_;

  private:
    static constexpr double NUM_LIMITS = 1.e-3;  ///< Numerical limits for sanity comparisons (MeV/mm-level)
    const Limits x_validity_range_{0., 1.};
    double kin_prefactor_{1.};
  };
}  // namespace cepgen::proc

#endif
