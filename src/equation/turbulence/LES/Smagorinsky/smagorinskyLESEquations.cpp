#include "smagorinskyLESEquations.h"
#include "mesh.h"
#include "messager.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace accel
{

smagorinskyLESEquations::smagorinskyLESEquations(realm* realm)
    : equation("Smagorinsky LES"), SmagorinskyLESModel(realm)
{
}

bool smagorinskyLESEquations::isConverged() const
{
    return true;
}

void smagorinskyLESEquations::setup()
{
    equation::isCreated_ = true;
}

void smagorinskyLESEquations::initialize()
{
    const char* perturbationText =
        std::getenv("OPENACCEL_CHANNEL_PERTURBATION");

    const scalar amplitude =
        perturbationText ? std::atof(perturbationText) : 0.0;

#if SPATIAL_DIM == 3
    if (amplitude > 0.0)
    {
        // Periodic-channel dimensions for the Re_tau = 180 case.
        constexpr scalar pi = 3.14159265358979323846;
        constexpr scalar Lx = 0.1256637061435917;
        constexpr scalar Ly = 0.02;
        constexpr scalar Lz = 0.04188790204786391;

        constexpr scalar kx = 2.0 * pi / Lx;
        constexpr scalar kz = 2.0 * pi / Lz;

        stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
        stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

        STKScalarField* USTKFieldPtr = this->URef().stkFieldPtr();

        const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, "coordinates");

        if (!USTKFieldPtr || !coordsSTKFieldPtr)
        {
            errorMsg(
                "Smagorinsky LES channel perturbation: velocity or "
                "coordinate field is unavailable");
        }

        for (const auto& domain : this->domainVector_)
        {
            stk::mesh::Selector selector =
                metaData.universal_part() &
                stk::mesh::selectUnion(
                    domain->zonePtr()->interiorParts());

            const stk::mesh::BucketVector& nodeBuckets =
                bulkData.get_buckets(
                    stk::topology::NODE_RANK, selector);

            for (auto ib = nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& bucket = **ib;

                for (stk::mesh::Bucket::size_type k = 0;
                     k < bucket.size();
                     ++k)
                {
                    const stk::mesh::Entity node = bucket[k];

                    const scalar* coordinates =
                        stk::mesh::field_data(
                            *coordsSTKFieldPtr, node);

                    scalar* velocity =
                        stk::mesh::field_data(
                            *USTKFieldPtr, node);

                    const scalar x = coordinates[0];
                    const scalar y = coordinates[1];
                    const scalar z = coordinates[2];

                    const scalar eta = y / Ly;
                    const scalar wallEnvelope =
                        std::sin(pi * eta) *
                        std::sin(pi * eta);
                    const scalar wallDerivativeShape =
                        std::sin(2.0 * pi * eta);

                    const scalar sinX = std::sin(kx * x);
                    const scalar cosX = std::cos(kx * x);
                    const scalar sinZ = std::sin(kz * z);
                    const scalar cosZ = std::cos(kz * z);

                    // Divergence-free, periodic perturbation.
                    const scalar uPrime =
                        amplitude * wallDerivativeShape * cosX;

                    const scalar vPrime =
                        2.0 * amplitude * Ly * wallEnvelope *
                        (sinX / Lx - sinZ / Lz);

                    const scalar wPrime =
                        -amplitude * wallDerivativeShape * cosZ;

                    velocity[0] += uPrime;
                    velocity[1] += vPrime;
                    velocity[2] += wPrime;
                }
            }
        }

        stk::mesh::communicate_field_data(
            bulkData, {USTKFieldPtr});

        if (messager::master())
        {
            std::cout
                << "[Smagorinsky LES] Applied periodic channel "
                << "velocity perturbation with amplitude "
                << amplitude << std::endl;
        }
    }
#else
    (void)amplitude;
#endif

    equation::isInitialized_ = true;
}

void smagorinskyLESEquations::postInitialize()
{
}

void smagorinskyLESEquations::preTimeStep()
{
}

void smagorinskyLESEquations::preSolve()
{
    FOREACH_DOMAIN(updateTurbulentDynamicViscosity);
    FOREACH_DOMAIN(updateEffectiveDynamicViscosity);
}

void smagorinskyLESEquations::solve()
{
}

} // namespace accel
