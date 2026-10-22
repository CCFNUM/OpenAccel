#pragma once

namespace accel
{

class realm;
class masking;

void declareImmersedDistance(realm& realm);
void computeImmersedDistance(realm& realm, const masking& masks);

}