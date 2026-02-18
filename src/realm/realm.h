#ifndef REALM_H
#define REALM_H

#include "types.h"

// fields
// issue: https://gitlab.cc-fnum.ch/cfd-code/cv-fem/accel-stk/-/issues/1
// resolving this issue would allow to just include something like this:
// #include "scalarNodeField.h"
// #include "vectorNodeField.h"
#include "compressibility.h"
#include "density.h"
#include "displacement.h"
#include "dynamicViscosity.h"
#include "heatFlowRate.h"
#include "massFlowRate.h"
#include "pressure.h"
#include "smRealm.h"
#include "specificEnthalpy.h"
#include "specificHeatCapacity.h"
#include "specificTotalEnthalpy.h"
#include "temperature.h"
#include "thermalConductivity.h"
#include "thermalExpansivity.h"
#include "turbRealm.h"
#include "velocity.h"
#include "volumeFraction.h"

namespace accel
{

class simulation;
class mesh;

class realm
{
private:
    // friend access
    friend class fieldBroker;

    simulation* simulationPtr_;

    mesh* meshPtr_;

    std::string name_;

    // Field declarations must be private

    std::unique_ptr<velocity> U_;

    std::unique_ptr<simpleVectorField> Ur_;

    std::unique_ptr<density> rho_;

    std::vector<std::unique_ptr<density>> rhoVector_; // for every phase

    std::unique_ptr<massFlowRate> mDot_;

    std::vector<std::unique_ptr<massFlowRate>> mDotVector_; // for every phase

    std::unique_ptr<pressure> p_;

    std::unique_ptr<simpleScalarField> p0_;

    std::unique_ptr<simpleScalarField> Ma_;

    std::unique_ptr<simpleScalarField> Co_;

    std::unique_ptr<dynamicViscosity> mu_;

    std::vector<std::unique_ptr<dynamicViscosity>> muVector_; // for every phase

    std::unique_ptr<temperature> T_;

    std::unique_ptr<simpleScalarField> T0_;

    std::unique_ptr<specificEnthalpy> h_;

    std::unique_ptr<specificTotalEnthalpy> h0_;

    std::unique_ptr<specificHeatCapacity> cp_;

    std::vector<std::unique_ptr<specificHeatCapacity>>
        cpVector_; // for every phase

    std::unique_ptr<thermalConductivity> lambda_;

    std::vector<std::unique_ptr<thermalConductivity>>
        lambdaVector_; // for every phase

    std::unique_ptr<thermalExpansivity> beta_;

    std::vector<std::unique_ptr<thermalExpansivity>>
        betaVector_; // for every phase

    std::unique_ptr<compressibility> psi_;

    std::vector<std::unique_ptr<compressibility>> psiVector_; // for every phase

    std::vector<std::unique_ptr<volumeFraction>>
        alphaVector_; // for every phase

    std::vector<std::unique_ptr<simpleVectorField>>
        nHatVector_; // for every phase

    std::vector<std::unique_ptr<simpleScalarField>>
        alphaSmoothVector_; // for every phase

    std::vector<std::unique_ptr<simpleScalarField>>
        rhsSmoothVector_; // for every phase

    std::unique_ptr<heatFlowRate> qDot_;

    std::unique_ptr<turbRealm> tRealm_;

    std::unique_ptr<nodeScalarField> yScale_;

    std::unique_ptr<simpleScalarField> minimumDistanceToWall_;

    // Mesh motion fields

    // total displacement from the original coordinates
    std::unique_ptr<nodeVectorField> Dt_ = nullptr;

    // displacement from the previous coordinates
    std::unique_ptr<displacement> D_ = nullptr;

    // mesh velocity
    std::unique_ptr<nodeVectorField> Um_ = nullptr;

    // divergence of mesh velocity
    std::unique_ptr<nodeScalarField> divUm_ = nullptr;

    std::unique_ptr<smRealm> smRealm_;

public:
    // field identifiers to be used with STK mesh field queries

    static constexpr char U_ID[] = "velocity";

    static constexpr char Ur_ID[] = "relative_velocity";

    static constexpr char rho_ID[] = "density";

    static constexpr char mDot_ID[] = "mass_flow_rate";

    static constexpr char p_ID[] = "pressure";

    static constexpr char mu_ID[] = "dynamic_viscosity";

    static constexpr char T_ID[] = "temperature";

    static constexpr char h_ID[] = "specific_enthalpy";

    static constexpr char h0_ID[] = "specific_total_enthalpy";

    static constexpr char cp_ID[] = "specific_heat_capacity";

    static constexpr char lambda_ID[] = "thermal_conductivity";

    static constexpr char beta_ID[] = "thermal_expansivity";

    static constexpr char psi_ID[] = "compressibility";

    static constexpr char alpha_ID[] = "volume_fraction";

    static constexpr char alpha_smooth_ID[] = "smoothed_volume_fraction";

    static constexpr char rhs_smooth_ID[] = "smoothed_rhs";

    static constexpr char qDot_ID[] = "heat_flow_rate";

    // other field identifiers

    static constexpr char p0_ID[] = "total_pressure";

    static constexpr char T0_ID[] = "total_temperature";

    static constexpr char Co_ID[] = "courant_number";

    static constexpr char Ma_ID[] = "mach_number";

    static constexpr char nHat_ID[] = "interface_normal";

    static constexpr char yScale_ID[] = "wall_scale";

    static constexpr char yMin_ID[] = "minimum_distance_to_wall";

    // Mesh motion identifiers

    static constexpr char Dt_ID[] = "displacement_mesh"; // total displacement

    static constexpr char D_ID[] = "displacement";

    static constexpr char Um_ID[] = "velocity_mesh";

    static constexpr char divUm_ID[] = "velocity_div_mesh";

    // Constructors

    realm(simulation* simulationPtr, std::string name);

    // Operations

    void initialize();

    void registerRestartField(const std::string& fieldName);

    // Access

    std::string name() const
    {
        return name_;
    };

    mesh* meshPtr();

    const mesh* meshPtr() const;

    mesh& meshRef();

    const mesh& meshRef() const;

    simulation* simulationPtr();

    const simulation* simulationPtr() const;

    simulation& simulationRef();

    const simulation& simulationRef() const;
};

} // namespace accel

#endif // REALM_H
