// File : turbRealm.h
// Created : Tue Jun 11 2024 15:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Turbulence field container for transport variables and wall
// functions
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBREALM_H
#define TURBREALM_H

#include "transitionOnsetReynoldsNumber.h"
#include "turbulentDissipationRate.h"
#include "turbulentEddyFrequency.h"
#include "turbulentIntermittency.h"
#include "turbulentKineticEnergy.h"
#include "turbulentThermalConductivity.h"
#include "turbulentViscosity.h"

namespace accel
{

class realm;
class fieldBroker;

class turbRealm
{
private:
    // friend access
    friend class realm;
    friend class fieldBroker;

    // Transport fields

    std::unique_ptr<turbulentKineticEnergy> k_;

    std::unique_ptr<turbulentEddyFrequency> omega_;

    std::unique_ptr<turbulentDissipationRate> epsilon_;

    std::unique_ptr<transitionOnsetReynoldsNumber> ReTheta_;

    std::unique_ptr<turbulentIntermittency> gamma_;

    // Turbulent properties

    std::unique_ptr<turbulentViscosity> mut_;

    std::unique_ptr<turbulentThermalConductivity> lambdat_;

    std::unique_ptr<nodeScalarField> muEff_;

    std::unique_ptr<nodeScalarField> lambdaEff_;

    std::unique_ptr<simpleScalarField> fOneBlending_;

    std::unique_ptr<simpleScalarField> tkeProduction_;

    // Wall-function related

    std::unique_ptr<simpleScalarField> uTau_;

    std::unique_ptr<simpleScalarField> yPlus_;

    std::unique_ptr<simpleScalarField> uPlus_;

    std::unique_ptr<simpleScalarField> yStar_;

    std::unique_ptr<simpleScalarField> uStar_;

    std::unique_ptr<simpleScalarField> duPlusdyPlus_;

    std::unique_ptr<simpleScalarField> uWallCoeffs_;

    std::unique_ptr<simpleVectorField> wallShearStress_;

    std::unique_ptr<simpleScalarField> TPlus_;

    std::unique_ptr<simpleScalarField> TWallCoeffs_;

public:
    // field identifiers to be used with STK mesh field queries

    static constexpr char k_ID[] = "turbulent_kinetic_energy";

    static constexpr char omega_ID[] = "turbulent_eddy_frequency";

    static constexpr char epsilon_ID[] = "turbulent_dissipation_rate";

    static constexpr char ReTheta_ID[] = "transition_onset_reynolds_number";

    static constexpr char gamma_ID[] = "turbulent_intermittency";

    static constexpr char mut_ID[] = "turbulent_viscosity";

    static constexpr char lambdat_ID[] = "turbulent_thermal_conductivity";

    static constexpr char muEff_ID[] = "effective_viscosity";

    static constexpr char lambdaEff_ID[] = "effective_thermal_conductivity";

    static constexpr char F1_ID[] = "f_one_blending";

    static constexpr char Pk_ID[] = "tke_production";

    static constexpr char yPlus_ID[] = "y_plus";

    static constexpr char uPlus_ID[] = "u_plus";

    static constexpr char yStar_ID[] = "y_star";

    static constexpr char uStar_ID[] = "u_star";

    static constexpr char uTau_ID[] = "wall_friction_velocity";

    static constexpr char duPlusdyPlus_ID[] = "du_plus_dy_plus";

    static constexpr char uWallCoeffs_ID[] = "u_wall_coeffs";

    static constexpr char wallShearStress_ID[] = "wall_shear_stress";

    static constexpr char TPlus_ID[] = "T_plus";

    static constexpr char TWallCoeffs_ID[] = "T_wall_coeffs";

    // Constructors

    turbRealm() = default;
};

} // namespace accel

#endif // TURBREALM_H
