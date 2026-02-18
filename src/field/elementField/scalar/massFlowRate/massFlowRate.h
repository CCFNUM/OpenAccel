// File : massFlowRate.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Element scalar field for mass flow rate with divergence and
// boundary fractions
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MASSFLOWRATE_H
#define MASSFLOWRATE_H

#include "elementScalarField.h"

namespace accel
{

template <size_t N, size_t M>
class nodeField;

class massFlowRate : public elementScalarField
{
protected:
    std::unique_ptr<nodeField<1, 1>> divFieldPtr_ = nullptr;

    std::vector<std::vector<scalar>> sideMassFlowRateFraction_;

public:
    // Constructors

    massFlowRate(realm* realmPtr,
                 const std::string name,
                 unsigned numberOfStates);

    // Methods

    void registerFractionSideField(label iZone, label iBoundary);

    // Access

    nodeField<1, 1>& divRef();

    const nodeField<1, 1>& divRef() const;

    // Getters and Setters

    void setF(label iZone, label iBoundary, scalar F)
    {
        sideMassFlowRateFraction_[iZone][iBoundary] = F;
    }

    scalar F(label iZone, label iBoundary) const
    {
        return sideMassFlowRateFraction_[iZone][iBoundary];
    }
};

} // namespace accel

#endif // MASSFLOWRATE_H
