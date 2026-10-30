#pragma once

namespace accel
{

class realm;
class masking;

void declareImmersedDistance(realm& realm);
struct ImmersedDistanceTimings
{
    double distance = 0.0;
    double sign = 0.0;
    double combined2D = 0.0;
    double field = 0.0;
};

void computeImmersedDistance(realm& realm, const masking& masks,
                             ImmersedDistanceTimings* timings = nullptr);

}