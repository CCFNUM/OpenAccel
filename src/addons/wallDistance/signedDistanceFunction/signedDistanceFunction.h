// File       : signedDistanceFunction.h
// Description: Wall distance computation using the SSDF library

#ifndef SIGNEDDISTANCEFUNCTION_H
#define SIGNEDDISTANCEFUNCTION_H

// code
#include "types.h"

namespace accel
{

class wallDistance;

class signedDistanceFunction
{
private:
    wallDistance* wallDistancePtr_;

    void computeWallDistance_();

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
