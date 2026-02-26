// File       : smRealm.h
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Solid mechanics field container for material properties and
// stress-strain
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SMREALM_H
#define SMREALM_H

// code
#include "poissonRatio.h"
#include "youngModulus.h"

namespace accel
{

class realm;
class fieldBroker;

class smRealm
{
private:
    // friend access
    friend class realm;
    friend class fieldBroker;

    // solid mechanics properties

    std::unique_ptr<youngModulus> E_;

    std::unique_ptr<poissonRatio> nu_;

    // solid mechanics post-processing fields

    // Cauchy stress tensor
    std::unique_ptr<simpleTensorField> sigma_ = nullptr;

    // strain tensor
    std::unique_ptr<simpleTensorField> epsilon_ = nullptr;

public:
    // field identifiers to be used with STK mesh field queries
    static constexpr char E_ID[] = "young_modulus";

    static constexpr char nu_ID[] = "poisson_ratio";

    static constexpr char sigma_ID[] = "stress";

    static constexpr char epsilon_ID[] = "strain";

    // Constructors

    smRealm() = default;
};

} // namespace accel

#endif // SMREALM_H
