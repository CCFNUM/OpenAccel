// File       : ggiInterfaceSideInfo.h
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: GGI (General Grid Interface) rasterization-based implementation
//              of interfaceSideInfo. One ggiInfo = one CS connection (one
//              current IP + one opposing face + area fraction).
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef GGIINTERFACESIDEINFO_H
#define GGIINTERFACESIDEINFO_H

#ifdef HAS_INTERFACE

// code
#include "ggiInfo.h"
#include "interfaceSideInfo.h"

namespace accel
{

class ggiInterfaceSideInfo : public interfaceSideInfo
{
public:
    ggiInterfaceSideInfo(interface* interfPtr,
                         bool isMasterSide,
                         const stk::mesh::PartVector currentPartVec,
                         const stk::mesh::PartVector opposingPartVec,
                         interfaceModelOption option,
                         label bitmapResolution,
                         scalar discernibleFraction,
                         bool verbose,
                         bool writeDiagnosisFile,
                         const std::string name);

    ~ggiInterfaceSideInfo();

    // -----------------------------------------------------------------------
    // virtual overrides
    // -----------------------------------------------------------------------

    void setup() override;
    void initialize() override;
    void update() override;
    void completeSearch() override;
    void determineElemsToGhost() override;
    void provideDiagnosis() override;
    size_t errorCheck() override;

private:
    // -----------------------------------------------------------------------
    // Per-face collected data: filled in collectFaceData_(), used by
    // constructBoundingPoints_() and completeSearch().
    // -----------------------------------------------------------------------
    struct FaceData_
    {
        stk::mesh::Entity face;
        stk::mesh::Entity element;
        label faceOrdinal;
        uint64_t globalFaceId;
        stk::topology faceTopo;
        stk::topology elemTopo;
        MasterElement* meFC;
        MasterElement* meSCS;
        label numIps;
        label numNodes;
        std::vector<scalar> nodeCoords; // [numNodes × SPATIAL_DIM]
        std::vector<scalar> ipCoords;   // [numIps   × SPATIAL_DIM]
        // SCS polygon vertices for each IP: scsVerts[ip][vertex 0-3][coord 0-2]
        // Vertex ordering (node, fwd_edge, centroid, bwd_edge):
        //   [0] = face node ip
        //   [1] = midpoint(node_ip, node_{ip+1})  (forward edge)
        //   [2] = face centroid
        //   [3] = midpoint(node_{ip-1}, node_ip)  (backward edge)
        std::vector<std::array<std::array<scalar, 3>, 4>> scsVerts;
    };

    // -----------------------------------------------------------------------
    // Data members (GGI-specific only; shared members live on
    // interfaceSideInfo)
    // -----------------------------------------------------------------------

    label bitmapResolution_; // NBIT: pixel-grid dimension (default 128)
    scalar
        discernibleFraction_; // minimum area fraction to retain (default 0.01)
    bool verbose_; // if true, print stats to stdout after completeSearch
    bool writeDiagnosisFile_; // if true, write per-rank GGI diagnosis files

    // Collected per-face data (filled by collectFaceData_())
    std::vector<FaceData_> faceDataVec_;

    // -----------------------------------------------------------------------
    // Private methods
    // -----------------------------------------------------------------------

    void deleteGgiInfo_();
    void reset_();
    void collectFaceData_();         // collect per-face data into faceDataVec_
    void constructBoundingPoints_(); // one circumsphere per current face
    void constructBoundingBoxes_();  // one AABB per opposing face
    void doSearch_();
};

} // namespace accel

#endif /* HAS_INTERFACE */
#endif
