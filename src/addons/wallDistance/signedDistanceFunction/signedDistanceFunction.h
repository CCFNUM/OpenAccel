// File       : signedDistanceFunction.h
// Description: Wall distance computation using the SSDF library

#ifndef SIGNEDDISTANCEFUNCTION_H
#define SIGNEDDISTANCEFUNCTION_H

// code
#include "masking.h"
#include "types.h"

#include <string>
#include <vector>

namespace accel
{

class wallDistance;

class signedDistanceFunction
{
private:
    wallDistance* wallDistancePtr_;

    void computeWallDistance_();
    void computeMaskedImmersedDistance_();

    static constexpr char immersedDistanceFieldName_[] =
        "signed_distance";

    void declareImmersedDistanceField_();

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
