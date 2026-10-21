// File       : multiphaseModel.h
// Created    : Sun Jan 26 2025 22:02:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Base multiphase model managing phase definitions and volume
//              fractions
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MULTIPHASEMODEL_H
#define MULTIPHASEMODEL_H

// code
#include "flowModel.h"

namespace accel
{

struct phase
{
    phase() = default;

    phase(label index, std::string name) : index_(index), name_(name)
    {
    }

    phase(label index, std::string name, bool primaryPhase)
        : index_(index), name_(name), primaryPhase_(primaryPhase)
    {
    }

    // global index in global set of materials in the simulation
    label index_ = -1;

    std::string name_;

    bool primaryPhase_ = true;
};

class multiphaseModel : public flowModel
{
protected:
    std::vector<phase> phases_;

public:
    // Constructors

    multiphaseModel(realm* realm);

    // Enable public use

    using fieldBroker::alphaRef;

    // Methods

    label phaseIndex(label iPhase) const
    {
        assert(iPhase < nPhases());
        return phases_[iPhase].index_;
    }

    label nPhases() const
    {
        return phases_.size();
    }

    phase& phaseRef(label iPhase)
    {
        return phases_[iPhase];
    }

    const phase& phaseRef(label iPhase) const
    {
        return phases_[iPhase];
    }

    // Interphase mass transfer (fields per fluid pair, liquid -> vapor > 0)

    bool hasMassTransfer() const
    {
        return !mdotSTKFieldPtrs_.empty();
    }

    bool hasMassTransfer(const domain* domain) const
    {
        for (const auto& fpm : domain->fluidPairModels())
        {
            if (fpm.massTransfer_.option_ != massTransferModelOption::none)
            {
                return true;
            }
        }
        return false;
    }

    const STKScalarField*
    massTransferRateSTKFieldPtr(const fluidPairModel& fpm) const
    {
        return mdotSTKFieldPtrs_.at(
            {fpm.massTransfer_.liquidIndex_, fpm.massTransfer_.vaporIndex_});
    }

    const STKScalarField*
    massTransferRateAlphaCoeffSTKFieldPtr(const fluidPairModel& fpm) const
    {
        return dmdotdalphaSTKFieldPtrs_.at(
            {fpm.massTransfer_.liquidIndex_, fpm.massTransfer_.vaporIndex_});
    }

    const STKScalarField*
    massTransferRatePressureCoeffSTKFieldPtr(const fluidPairModel& fpm) const
    {
        return dmdotdpSTKFieldPtrs_.at(
            {fpm.massTransfer_.liquidIndex_, fpm.massTransfer_.vaporIndex_});
    }

    // restore the rate on restart (zero otherwise)
    void initializeMassTransferRate(const std::shared_ptr<domain> domain);

    // raw rate from p, alpha, rho, under-relaxed against the stored one
    void updateMassTransferRate(const std::shared_ptr<domain> domain);

    // follow the last pressure change with the linearized rate
    void correctMassTransferRate(const std::shared_ptr<domain> domain);

protected:
    // mass transfer rate and its sensitivities d/d alpha_v, d/d p, keyed by
    // the global material indices of the liquid and vapor of the pair
    std::map<std::pair<label, label>, STKScalarField*> mdotSTKFieldPtrs_;
    std::map<std::pair<label, label>, STKScalarField*> dmdotdalphaSTKFieldPtrs_;
    std::map<std::pair<label, label>, STKScalarField*> dmdotdpSTKFieldPtrs_;

    // requires the phases to be collected
    void setupMassTransfer_(realm* realm);
};

} /* namespace accel */

#endif // MULTIPHASEMODEL_H
