// File       : solidMechanicsModel.h
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Solid mechanics model for displacement, stress, and strain
// computations
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SOLIDMECHANICSMODEL_H
#define SOLIDMECHANICSMODEL_H

// code
#include "model.h"

namespace accel
{

class solidMechanicsModel : public model
{
private:
    // update

    void updateDisplacementSideFields_(const std::shared_ptr<domain> domain);

    void updateDisplacementBoundarySideFieldTraction_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

#ifdef HAS_INTERFACE
    void updateDisplacementInterfaceSideFieldTraction_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr);
#endif /* HAS_INTERFACE */

public:
    solidMechanicsModel(realm* realm);

    using fieldBroker::DRef;
    using fieldBroker::ERef;
    using fieldBroker::nuRef;
    using fieldBroker::pRef;
    using fieldBroker::strainRef;
    using fieldBroker::stressRef;
    using fieldBroker::wallShearStressRef;

    // setup

    void setupDisplacement(const std::shared_ptr<domain> domain) override;

    // update

    void updateDisplacement(const std::shared_ptr<domain> domain) override;

protected:
    void updateStressAndStrain_(const std::shared_ptr<domain> domain);

    std::string
    getCoordinatesID_(const std::shared_ptr<domain> domain) const override
    {
        bool useOrig = false;
        if (domain->solidMechanics_.formulation_ ==
                kinematicFormulationType::totalLagrangian &&
            this->controlsRef().isTransient())
        {
            useOrig = true;
        }
        return useOrig ? mesh::original_coordinates_ID : mesh::coordinates_ID;
    }

    std::string
    getDualNodalVolumeID_(const std::shared_ptr<domain> domain) const override
    {
        bool useOrig = false;
        if (domain->solidMechanics_.formulation_ ==
                kinematicFormulationType::totalLagrangian &&
            this->controlsRef().isTransient())
        {
            useOrig = true;
        }
        return useOrig ? mesh::original_dual_nodal_volume_ID
                       : mesh::dual_nodal_volume_ID;
    }

    std::string
    getExposedAreaVectorID_(const std::shared_ptr<domain> domain) const override
    {
        bool useOrig = false;
        if (domain->solidMechanics_.formulation_ ==
                kinematicFormulationType::totalLagrangian &&
            this->controlsRef().isTransient())
        {
            useOrig = true;
        }
        return useOrig ? mesh::original_exposed_area_vector_ID
                       : mesh::exposed_area_vector_ID;
    }
};

} /* namespace accel */

#endif // SOLIDMECHANICSMODEL_H
