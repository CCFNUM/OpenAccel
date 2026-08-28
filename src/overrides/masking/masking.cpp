// File       : masking.cpp
// Created    : Tue Aug 11 2026
// Author     : Mhamad Mahdi Alloush
// Description: Manager of the solid bodies immersed in the fluid mesh
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "masking.h"
#include "domain.h"
#include "macros.h"
#include "messager.h"
#include "realm.h"
#include "simulation.h"

#include <graph/layers.hpp>

namespace accel
{

masking::masking(realm* realmPtr) : fieldBroker(realmPtr)
{
}

void masking::read(const YAML::Node& maskingArrayNode)
{
    label id = 1;
    for (const auto& regionNode : maskingArrayNode)
    {
        auto regionPtr = std::make_unique<maskedRegion>(realmPtr_, id++);
        regionPtr->read(regionNode);
        regionVector_.push_back(std::move(regionPtr));
    }

    if (regionVector_.empty())
    {
        errorMsg("maskedRegion: no body defined");
    }
}

void masking::setup()
{
    auto& metaData = meshRef().metaDataRef();

    const auto selInterior =
        stk::mesh::selectUnion(meshRef().interiorActiveParts());

    for (const char* name : {layer_ID, id_ID, distance_ID})
    {
        auto& stkFieldRef =
            metaData.declare_field<scalar>(stk::topology::NODE_RANK, name);
        stk::mesh::put_field_on_mesh(stkFieldRef, selInterior, 1, nullptr);
        stk::io::set_field_output_type(stkFieldRef,
                                       stk::io::FieldOutputType::SCALAR);
    }
}

void masking::initialize()
{
    for (auto& regionPtr : regionVector_)
    {
        regionPtr->initialize();
    }

    buildGraph_();
    classify_();
}

void masking::update()
{
    buildGraph_();
    classify_();
}

void masking::classify_()
{
    const label npoints = graphNodeCount();

    std::vector<uint8_t> seeds(npoints, 255);
    std::vector<uint8_t> layers(npoints, 255);

    // the layer a node ends up in and the region that put it there
    std::vector<uint8_t> bestLayer(npoints, 255);
    std::vector<label> bestId(npoints, 0);

    for (label r = 0; r < regionCount(); ++r)
    {
        maskedRegion& region = regionRef(r);

        region.markCovered(seeds);

        // one layer past the reference layer is all the wall treatment needs
        computeLayers_(seeds, layers, region.referenceLayer() + 1);

        region.takeLayers(layers);

        for (label i = 0; i < npoints; ++i)
        {
            if (layers[i] < bestLayer[i])
            {
                bestLayer[i] = layers[i];
                bestId[i] = region.id();
            }
        }
    }

    STKScalarField* layerField = meshRef().metaDataRef().get_field<scalar>(
        stk::topology::NODE_RANK, layer_ID);
    STKScalarField* idField = meshRef().metaDataRef().get_field<scalar>(
        stk::topology::NODE_RANK, id_ID);

    for (label i = 0; i < npoints; ++i)
    {
        *stk::mesh::field_data(*layerField, graphNodes_[i]) =
            static_cast<scalar>(bestLayer[i]);
        *stk::mesh::field_data(*idField, graphNodes_[i]) =
            static_cast<scalar>(bestId[i]);
    }

    collect_();
}

label masking::graphIndexOf(const stk::mesh::Entity node) const
{
    const auto local = meshRef().bulkDataRef().local_id(node);

    return (local < graphIndex_.size()) ? graphIndex_[local] : -1;
}

void masking::buildGraph_()
{
    const auto& bulkData = meshRef().bulkDataRef();
    const auto& metaData = meshRef().metaDataRef();

    // the graph spans every masked domain, so it is shared by all regions
    stk::mesh::PartVector parts;
    for (label i = 0; i < regionCount(); ++i)
    {
        for (auto part : regionRef(i).domainParts())
        {
            parts.push_back(part);
        }
    }
    const stk::mesh::Selector selDomains = stk::mesh::selectUnion(parts);

    graphNodes_.clear();
    graphIndex_.assign(bulkData.get_size_of_entity_index_space(), -1);

    for (const stk::mesh::Bucket* bucketPtr :
         bulkData.get_buckets(stk::topology::NODE_RANK, selDomains))
    {
        for (const stk::mesh::Entity node : *bucketPtr)
        {
            graphIndex_[bulkData.local_id(node)] =
                static_cast<label>(graphNodes_.size());
            graphNodes_.push_back(node);
        }
    }

    // ssdf::layers takes a single node count per element: a shorter element is
    // padded by repeating its last node, which leaves its neighbours unchanged
    nodesPerElement_ = 0;
    label elementCount = 0;
    for (const stk::mesh::Bucket* bucketPtr :
         bulkData.get_buckets(stk::topology::ELEMENT_RANK, selDomains))
    {
        nodesPerElement_ = std::max<label>(
            nodesPerElement_,
            static_cast<label>(bucketPtr->topology().num_nodes()));
        elementCount += static_cast<label>(bucketPtr->size());
    }

    elementNodes_.assign(nodesPerElement_, std::vector<label>(elementCount, 0));

    label e = 0;
    for (const stk::mesh::Bucket* bucketPtr :
         bulkData.get_buckets(stk::topology::ELEMENT_RANK, selDomains))
    {
        for (const stk::mesh::Entity element : *bucketPtr)
        {
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
            const unsigned numNodes = bulkData.num_nodes(element);

            for (label i = 0; i < nodesPerElement_; ++i)
            {
                const unsigned n = std::min<unsigned>(i, numNodes - 1);
                elementNodes_[i][e] = graphIndexOf(nodes[n]);
            }
            ++e;
        }
    }

    elementNodeColumns_.resize(nodesPerElement_);
    for (label i = 0; i < nodesPerElement_; ++i)
    {
        elementNodeColumns_[i] = elementNodes_[i].data();
    }
}

void masking::computeLayers_(const std::vector<uint8_t>& seeds,
                             std::vector<uint8_t>& layers,
                             label maxLayers)
{
    const auto& bulkData = meshRef().bulkDataRef();
    const label npoints = graphNodeCount();
    const label nelements =
        nodesPerElement_ > 0 ? static_cast<label>(elementNodes_[0].size()) : 0;

    layers = seeds;

    const bool hasSeed =
        std::find(layers.begin(), layers.end(), 0) != layers.end();

    if (hasSeed && nelements > 0)
    {
        ssdf::layers<label>(nelements,
                            static_cast<int>(nodesPerElement_),
                            elementNodeColumns_.data(),
                            npoints,
                            layers.data(),
                            static_cast<int>(maxLayers));
    }
    else
    {
        std::fill(layers.begin(), layers.end(), 255);
    }

    if (!messager::parallel() || nelements == 0)
    {
        return;
    }

    // a region may be cut by a rank boundary: exchange the layer of the shared
    // nodes and continue the walk from whatever arrived, until nothing moves
    STKScalarField* layerField = meshRef().metaDataRef().get_field<scalar>(
        stk::topology::NODE_RANK, layer_ID);

    ptrdiff_t* n2ePtr = nullptr;
    label* elementIndex = nullptr;
    ssdf::create_n2e<label, ptrdiff_t, label>(
        nelements,
        npoints,
        static_cast<int>(nodesPerElement_),
        elementNodeColumns_.data(),
        &n2ePtr,
        &elementIndex);

    std::vector<label> queue(npoints);

    for (label sweep = 0; sweep < messager::nProcs() + 1; ++sweep)
    {
        for (label i = 0; i < npoints; ++i)
        {
            *stk::mesh::field_data(*layerField, graphNodes_[i]) =
                static_cast<scalar>(layers[i]);
        }

        stk::mesh::parallel_min(bulkData, {layerField});

        int queueSize = 0;
        for (label i = 0; i < npoints; ++i)
        {
            const uint8_t received = static_cast<uint8_t>(
                *stk::mesh::field_data(*layerField, graphNodes_[i]));

            if (received < layers[i])
            {
                layers[i] = received;
                queue[queueSize++] = i;
            }
        }

        int changed = queueSize;
        MPI_Allreduce(
            MPI_IN_PLACE, &changed, 1, MPI_INT, MPI_SUM, messager::comm());

        if (changed == 0)
        {
            break;
        }

        int cursor = 0;
        ssdf::layers_iterative<label>(static_cast<int>(nodesPerElement_),
                                      elementNodeColumns_.data(),
                                      layers.data(),
                                      static_cast<int>(maxLayers),
                                      n2ePtr,
                                      elementIndex,
                                      &cursor,
                                      &queueSize,
                                      static_cast<int>(queue.size()),
                                      queue.data());
    }

    free(n2ePtr);
    free(elementIndex);
}

void masking::collect_()
{
    coveredNodes_.clear();
    ringNodes_.clear();
    sx_.clear();
    sy_.clear();
    sz_.clear();
    s0_.clear();
    s1_.clear();
    s2_.clear();

    for (const auto& regionPtr : regionVector_)
    {
        coveredNodes_.insert(coveredNodes_.end(),
                             regionPtr->coveredNodes().begin(),
                             regionPtr->coveredNodes().end());
        ringNodes_.insert(ringNodes_.end(),
                          regionPtr->ringNodes().begin(),
                          regionPtr->ringNodes().end());

        const label offset = static_cast<label>(sx_.size());

        sx_.insert(sx_.end(),
                   regionPtr->surfaceX().begin(),
                   regionPtr->surfaceX().end());
        sy_.insert(sy_.end(),
                   regionPtr->surfaceY().begin(),
                   regionPtr->surfaceY().end());
        sz_.insert(sz_.end(),
                   regionPtr->surfaceZ().begin(),
                   regionPtr->surfaceZ().end());

        for (label f = 0; f < regionPtr->facetCount(); ++f)
        {
            s0_.push_back(offset + regionPtr->facetNode0()[f]);
            s1_.push_back(offset + regionPtr->facetNode1()[f]);
            s2_.push_back(offset + regionPtr->facetNode2()[f]);
        }
    }
}

void masking::setAtCovered(STKScalarField& field, scalar value)
{
    for (const auto& regionPtr : regionVector_)
    {
        regionPtr->setAtCovered(field, value);
    }
}

void masking::copyAtCovered(STKScalarField& target,
                            const STKScalarField& source)
{
    for (const auto& regionPtr : regionVector_)
    {
        regionPtr->copyAtCovered(target, source);
    }
}

void masking::reapplyVelocity()
{
    auto& U = URef().stkFieldRef();

    for (const auto& regionPtr : regionVector_)
    {
        regionPtr->reapplyVelocity(U);
    }
}

void masking::reportForces() const
{
    if (!messager::master())
    {
        return;
    }

    for (const auto& regionPtr : regionVector_)
    {
        const scalar* f = regionPtr->force();
        std::cout << "Masked region `" << regionPtr->name() << "`: force = ["
                  << f[0] << ", " << f[1];
#if SPATIAL_DIM == 3
        std::cout << ", " << f[2];
#endif
        std::cout << "] [N]" << std::endl;
    }
}

} // namespace accel
