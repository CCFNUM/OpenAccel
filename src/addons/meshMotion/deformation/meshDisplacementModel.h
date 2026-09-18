// File       : meshDisplacementModel.h
// Created    : Mon Sep 14 2026
// Author     : Mhamad Mahdi Alloush
// Description: Mesh displacement field with its wall and interface motion,
//              shared by the displacement diffusion and Lithe backends
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MESHDISPLACEMENTMODEL_H
#define MESHDISPLACEMENTMODEL_H

// code
#include "model.h"

namespace accel
{

// wall/interface motion of the mesh displacement and its Dirichlet parts
class meshDisplacementModel : public model
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

    void updateDisplacementInterfaceSideFieldDeformation_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr);

    // Read the side's mesh_motion displacement data into its dataHandler
    void setupInterfaceSideMeshMotion_(interfaceSideInfo* interfaceSideInfoPtr,
                                       const YAML::Node& sideNode);

    // Impose the side's mesh_motion option (stationary/specified/periodic)
    void updateDisplacementInterfaceSidePrescribed_(
        const std::shared_ptr<domain> domain,
        interfaceSideInfo* interfaceSideInfoPtr);

    // Calculate force and moment over a specified patch
    void calculateSurfaceForceAndMoment_(const boundary* boundary,
                                         const utils::vector& center,
                                         utils::vector& force,
                                         utils::vector& moment);

protected:
    // parts whose displacement is prescribed (walls, moving sides, FSI sides)
    stk::mesh::PartVector collectDirichletBoundaryParts_(const domain* domain);

public:
    meshDisplacementModel(realm* realm);

    using fieldBroker::DRef;
    using fieldBroker::DtRef;
    using fieldBroker::pRef;
    using fieldBroker::wallShearStressRef;
    using fieldBroker::yMinRef;

    // setup

    void setupDisplacement(const std::shared_ptr<domain> domain) override;

    // initialize

    void initializeDisplacement(const std::shared_ptr<domain> domain) override;

    // update

    void updateDisplacement(const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // MESHDISPLACEMENTMODEL_H
