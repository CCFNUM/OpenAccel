// File       : displacementDiffusionModel.h
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Model for displacement boundary conditions and surface loads
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DISPLACEMENTDIFFUSIONMODEL_H
#define DISPLACEMENTDIFFUSIONMODEL_H

// code
#include "model.h"

namespace accel
{

class displacementDiffusionModel : public model
{
private:
    void updateDisplacementSideFields_(const std::shared_ptr<domain> domain);

    void updateDisplacementBoundarySideFieldSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateDisplacementBoundarySideFieldPeriodicDisplacement_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateDisplacementBoundarySideFieldRigidBodySolution_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

#ifdef HAS_INTERFACE
    void updateDisplacementInterfaceSideFieldDeformation_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr);
#endif /* HAS_INTERFACE */

    // Calculate force and moment over a specified patch
    void calculateSurfaceForceAndMoment_(const boundary* boundary,
                                         const utils::vector& center,
                                         utils::vector& force,
                                         utils::vector& moment);

public:
    displacementDiffusionModel(realm* realm);

    using fieldBroker::DRef;
    using fieldBroker::DtRef;
    using fieldBroker::pRef;
    using fieldBroker::wallShearStressRef;
    using fieldBroker::yMinRef;

    // initialize

    void initializeDisplacement(const std::shared_ptr<domain> domain) override;

    // update

    void updateDisplacement(const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // DISPLACEMENTDIFFUSIONMODEL_H
