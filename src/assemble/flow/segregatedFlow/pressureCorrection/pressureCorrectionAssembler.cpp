// File       : pressureCorrectionAssembler.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "pressureCorrectionAssembler.h"
#include "flowModel.h"

namespace accel
{

pressureCorrectionAssembler::pressureCorrectionAssembler(flowModel* model)
    : phiAssembler<1>(static_cast<fieldBroker*>(model)), model_(model)
{
}

void pressureCorrectionAssembler::adjustMatrixForPressureReference(
    const domain* domain,
    Context* ctx)
{
    assert(Base::BLOCKSIZE == 1);

    // processor rank of the zone partition carrying the reference node
    label refRank = domain->associatedPartitionRankForPressureLevelNode();

    if (refRank == messager::myProcNo())
    {
        // apply to linear system
        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        // reference node
        stk::mesh::EntityId gid = domain->pressureLevelNodeId();
        stk::mesh::Entity refNode = model_->meshRef().bulkDataRef().get_entity(
            stk::topology::NODE_RANK, gid);
        const label& lid = model_->meshRef().bulkDataRef().local_id(refNode);

        // get reference to diagonal and store it
        scalar& dia = A.dofDiag(lid);
        scalar d = dia;

        // zero-out the whole row of the pressure-correction equation at the
        // given node
        auto cols = A.rowCols(lid);
        auto vals = A.rowVals(lid);

        for (label jdx = 0; jdx < cols.size(); jdx++)
        {
            vals[jdx] = 0.0;
        }

        // restore diagonal
        dia = d;

        // grab the rhs of the pressure-correction equation
        scalar& rhs = b[lid];

        // set rhs
        rhs = 0;
    }
}

} // namespace accel
