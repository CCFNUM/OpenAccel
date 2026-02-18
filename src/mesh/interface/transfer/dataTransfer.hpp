// File       : dataTransfer.hpp
// Created    : Fri Nov 21 2025 14:01:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <class FROM, class TO>
void nonconformalDataTransfer::linearInterpolation<FROM, TO>::filter_to_nearest(
    EntityKeyMap& RangeToDomain,
    const MeshA& fromElem,
    MeshB& toPoints)
{
    const stk::mesh::BulkData& fromBulkData = fromElem.fromBulkData_;
    stk::mesh::BulkData& toBulkData = toPoints.toBulkData_;

    const stk::mesh::Field<scalar>* fromcoordinates = fromElem.fromcoordinates_;
    const stk::mesh::Field<scalar>* tocoordinates = toPoints.tocoordinates_;

    typedef typename EntityKeyMap::iterator iterator;
    typedef typename EntityKeyMap::const_iterator const_iterator;

    // some simple user diagnostics to let the user know if something bad is
    // happening
    scalar maxBestX = -std::numeric_limits<scalar>::max();
    size_t maxCandidateBoundingBox = 0;

    for (const_iterator current_key = RangeToDomain.begin();
         current_key != RangeToDomain.end();)
    {
        scalar bestX_ = std::numeric_limits<scalar>::max();

        const stk::mesh::EntityKey thePt = current_key->first;
        stk::mesh::Entity theNode = toBulkData.get_entity(thePt);
        // load nodal coordinates from node
        const scalar* tocoords = stk::mesh::field_data(*tocoordinates, theNode);

        std::pair<iterator, iterator> keys =
            RangeToDomain.equal_range(current_key->first);
        iterator nearest = keys.second;

        size_t candidateBoundingBoxSize = 0;
        for (iterator ii = keys.first; ii != keys.second; ++ii)
        {
            candidateBoundingBoxSize++;

            const stk::mesh::EntityKey theBox = ii->second;
            stk::mesh::Entity theElem = fromBulkData.get_entity(theBox);

            // extract master element from the bucket in which the element
            // resides
            const stk::mesh::Bucket& theBucket = fromBulkData.bucket(theElem);
            const stk::topology& theElemTopo = theBucket.topology();
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(theElemTopo);

            // load nodal coordinates from element
            stk::mesh::Entity const* elem_node_rels =
                fromBulkData.begin_nodes(theElem);
            const label num_nodes = fromBulkData.num_nodes(theElem);

            const label nodesPerElement = meSCS->nodesPerElement_;
            std::vector<scalar> theElementCoords(SPATIAL_DIM * nodesPerElement);

            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = elem_node_rels[ni];

                // load up vectors
                const scalar* fromcoords =
                    stk::mesh::field_data(*fromcoordinates, node);
                for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                {
                    const label offSet = j * nodesPerElement + ni;
                    theElementCoords[offSet] = fromcoords[j];
                }
            }

            std::vector<scalar> isoParCoords(SPATIAL_DIM);
            const scalar nearestDistance = meSCS->isInElement(
                &theElementCoords[0], &(tocoords[0]), &(isoParCoords[0]));
            if (nearestDistance < bestX_)
            {
                bestX_ = nearestDistance;
                toPoints.TransferInfo_[thePt] = isoParCoords;
                nearest = ii;
                maxBestX = std::max(maxBestX, bestX_);
            }
        }
        maxCandidateBoundingBox =
            std::max(maxCandidateBoundingBox, candidateBoundingBoxSize);

        current_key = keys.second;
        if (nearest != keys.first)
            RangeToDomain.erase(keys.first, nearest);
        if (nearest != keys.second)
            RangeToDomain.erase(++nearest, keys.second);
    }

    // parallel sum and output diagnostics
    scalar g_maxBestX = 0.0;
    size_t g_maxCandidateBoundngBox = 0;
    stk::all_reduce_max(fromElem.comm(), &maxBestX, &g_maxBestX, 1);
    stk::all_reduce_max(fromElem.comm(),
                        &maxCandidateBoundingBox,
                        &g_maxCandidateBoundngBox,
                        1);

    if (toBulkData.parallel_rank() == 0)
    {
        std::cout << std::endl;
        std::cout << "XFER::LinInterp::fine_search() Overview:" << std::endl;
        std::cout << "  Maximum normalized distance found is: " << g_maxBestX
                  << " (should be unity or less)" << std::endl;
        std::cout << "  Maximum number of candidate bounding boxes found for a "
                     "single point is: "
                  << g_maxCandidateBoundngBox << std::endl;
        std::cout << "  Should max normalized distance and/or candidate "
                     "bounding box size be too large, please check setup"
                  << std::endl;
    }
}

template <class FROM, class TO>
void nonconformalDataTransfer::linearInterpolation<FROM, TO>::apply(
    MeshB& ToPoints,
    const MeshA& FromElem,
    const EntityKeyMap& RangeToDomain)
{
    const stk::mesh::BulkData& fromBulkData = FromElem.fromBulkData_;
    stk::mesh::BulkData& toBulkData = ToPoints.toBulkData_;

    typename EntityKeyMap::const_iterator ii;
    for (ii = RangeToDomain.begin(); ii != RangeToDomain.end(); ++ii)
    {
        const stk::mesh::EntityKey thePt = ii->first;
        const stk::mesh::EntityKey theBox = ii->second;

        if (1 != ToPoints.TransferInfo_.count(thePt))
        {
            if (0 == ToPoints.TransferInfo_.count(thePt))
                throw std::runtime_error("Key not found in database");
            else
                throw std::runtime_error("Too many Keys found in database");
        }
        const std::vector<scalar>& isoParCoords_ =
            ToPoints.TransferInfo_[thePt];
        stk::mesh::Entity theNode = toBulkData.get_entity(thePt);
        stk::mesh::Entity theElem = fromBulkData.get_entity(theBox);

        const stk::mesh::Bucket& theBucket = fromBulkData.bucket(theElem);
        const stk::topology& theElemTopo = theBucket.topology();
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);

        stk::mesh::Entity const* elem_node_rels =
            fromBulkData.begin_nodes(theElem);
        const label num_nodes = fromBulkData.num_nodes(theElem);
        const label nodesPerElement = meSCS->nodesPerElement_;

        for (unsigned n = 0; n != FromElem.fromFieldVec_.size(); ++n)
        {
            // extract field
            const stk::mesh::FieldBase* toFieldBaseField =
                ToPoints.toFieldVec_[n];

            // extract field name
            const std::string fieldName = toFieldBaseField->name();

            // find any clipping
            scalar clipMin = std::numeric_limits<scalar>::lowest();
            scalar clipMax = std::numeric_limits<scalar>::max();
            std::map<std::string, std::pair<scalar, scalar>>::iterator itc =
                ToPoints.clipMap_.find(fieldName);
            if (itc != ToPoints.clipMap_.end())
            {
                clipMin = (*itc).second.first;
                clipMax = (*itc).second.second;
            }

            // FixMe: integers are problematic for now...
            const size_t sizeOfField =
                field_bytes_per_entity(*toFieldBaseField, theNode) /
                sizeof(scalar);
            std::vector<scalar> Coeff(nodesPerElement * sizeOfField);

            // now load the elemental values for future interpolation; fill in
            // connected nodes
            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = elem_node_rels[ni];

                const stk::mesh::FieldBase* fromFieldBaseField =
                    FromElem.fromFieldVec_[n];
                const scalar* theField =
                    (scalar*)stk::mesh::field_data(*fromFieldBaseField, node);

                for (size_t j = 0; j < sizeOfField; ++j)
                {
                    const label offSet = j * nodesPerElement + ni;
                    Coeff[offSet] = theField[j];
                }
            }

            scalar* toField =
                (scalar*)stk::mesh::field_data(*toFieldBaseField, theNode);
            if (!toField)
                throw std::runtime_error(
                    "Receiving field undefined on mesh object.");
            meSCS->interpolatePoint(
                sizeOfField, &isoParCoords_[0], &Coeff[0], toField);

            // clip it
            for (size_t j = 0; j < sizeOfField; ++j)
            {
                toField[j] = std::min(clipMax, std::max(toField[j], clipMin));
            }
        }
    }
}

} // namespace accel
