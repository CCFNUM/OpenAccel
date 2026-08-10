// File       : model.cpp
// Created    : Wed Jan 07 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "model.h"
#include "domain.h"
#include "messager.h"
#include "realm.h"

namespace accel
{

model::model(realm* realm) : fieldBroker(realm)
{
}

bool model::boundaryInputAreaAverage_(const std::shared_ptr<domain> domain,
                                      label iBoundary,
                                      inputData<1>& data,
                                      scalar& average,
                                      scalar& area)
{
    average = 0.0;
    area = 0.0;

    const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

    const auto& mesh = this->meshRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& areaField = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    const bool evaluateExpression = data.type() == inputDataType::expression;

    // exprtk context; z stays zero in 2D to avoid an unset variable
    typedef exprtk::symbol_table<scalar> symbol_table_t;
    typedef exprtk::expression<scalar> expression_t;
    typedef exprtk::parser<scalar> parser_t;

    symbol_table_t symbolTable;
    scalar t = this->controlsRef().time;
    scalar x = 0.0, y = 0.0, z = 0.0;
    expression_t expression;

    if (evaluateExpression)
    {
        symbolTable.add_constants();
        symbolTable.add_variable("t", t);
        symbolTable.add_variable("x", x);
        symbolTable.add_variable("y", y);
        symbolTable.add_variable("z", z);
        expression.register_symbol_table(symbolTable);

        parser_t parser;
        if (!parser.compile(data.expression()[0], expression))
        {
            errorMsg("Error in the expression provided at boundary `" +
                     boundaryRef.name() + "`: " + data.expression()[0]);
        }
    }

    // workspace
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_face_shape_function;
    std::array<scalar, SPATIAL_DIM> coordBip{};

    scalar localArea = 0.0;
    scalar localSum = 0.0;

    const auto sideSelector = metaData.locally_owned_part() &
                              stk::mesh::selectUnion(boundaryRef.parts());
    const auto& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), sideSelector);
    for (const auto* bucketPtr : sideBuckets)
    {
        const auto& bucket = *bucketPtr;

        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket.topology());
        const label nodesPerSide = bucket.topology().num_nodes();
        const label numScsBip = meFC->numIntPoints_;

        ws_coordinates.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);
        meFC->shape_fcn(&ws_face_shape_function[0]);

        for (const auto side : bucket)
        {
            const scalar* areav = stk::mesh::field_data(areaField, side);

            if (evaluateExpression)
            {
                const stk::mesh::Entity* sideNodeRels =
                    bulkData.begin_nodes(side);
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    const scalar* coords =
                        stk::mesh::field_data(coordsField, sideNodeRels[ni]);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        ws_coordinates[ni * SPATIAL_DIM + j] = coords[j];
                    }
                }
            }

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                scalar areaMagnitudeSquared = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar component = areav[ip * SPATIAL_DIM + j];
                    areaMagnitudeSquared += component * component;
                }
                const scalar areaMagnitude = std::sqrt(areaMagnitudeSquared);
                localArea += areaMagnitude;

                if (!evaluateExpression)
                {
                    continue;
                }

                // interpolate the boundary integration point coordinate
                coordBip.fill(0.0);
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = ws_face_shape_function[offSetSF_face + ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        coordBip[j] += r * ws_coordinates[ic * SPATIAL_DIM + j];
                    }
                }

                x = coordBip[0];
                y = coordBip[1];
#if SPATIAL_DIM == 3
                z = coordBip[2];
#endif
                localSum += expression.value() * areaMagnitude;
            }
        }
    }

    scalar localReduction[2] = {localArea, localSum};
    scalar globalReduction[2] = {0.0, 0.0};
    stk::all_reduce_sum(
        bulkData.parallel(), localReduction, globalReduction, 2);

    area = globalReduction[0];

    switch (data.type())
    {
        case inputDataType::constant:
            average = data.value()[0];
            return true;

        case inputDataType::timeTable:
            average = data.interpolate(this->controlsRef().time)[0];
            return true;

        case inputDataType::expression:
            if (area <= SMALL)
            {
                return false;
            }
            average = globalReduction[1] / area;
            return true;

        case inputDataType::profileData:
            {
                // initial guess only: samples are not area weighted
                const auto& values = data.scatterValues();
                if (values.empty())
                {
                    return false;
                }
                average =
                    std::accumulate(values.begin(), values.end(), scalar(0)) /
                    static_cast<scalar>(values.size());
                return true;
            }

        default:
            return false;
    }
}

} /* namespace accel */
