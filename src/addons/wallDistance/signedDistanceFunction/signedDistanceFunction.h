// File       : signedDistanceFunction.h
// Description: Domain-wall distance computation using SSDF

#ifndef SIGNEDDISTANCEFUNCTION_H
#define SIGNEDDISTANCEFUNCTION_H

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
    signedDistanceFunction(wallDistance* wallDistancePtr);

    void setup();

    void initialize();

    void update();
};

} // namespace accel

#endif // SIGNEDDISTANCEFUNCTION_H