// File       : signedDistanceFunction.h
// Description: Wall distance computation using the SSDF library

#ifndef SIGNEDDISTANCEFUNCTION_H
#define SIGNEDDISTANCEFUNCTION_H

// code
#include "masking.h"
#include "types.h"

namespace accel
{

class wallDistance;

class signedDistanceFunction
{
private:
    wallDistance* wallDistancePtr_;

    void computeWallDistance_();

    void computeImmersedSolidDistance_(
        const masking& masks,
        const std::vector<stk::mesh::Entity>& queryNodes,
        const std::vector<scalar>& x,
        const std::vector<scalar>& y,
        const std::vector<scalar>& z);

public:
    // Constructors
    signedDistanceFunction(wallDistance* wallDistancePtr);

    // Methods

    void setup();

    void initialize();

    void update();
};

} // namespace accel

#endif // SIGNEDDISTANCEFUNCTION_H
