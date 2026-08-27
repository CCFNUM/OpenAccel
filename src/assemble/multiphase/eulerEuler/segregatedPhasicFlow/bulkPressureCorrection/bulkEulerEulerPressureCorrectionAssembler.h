// File       : bulkEulerEulerPressureCorrectionAssembler.h
// Description: Common pressure correction assembled from all phase fluxes

#ifndef BULKEULEREULERPRESSURECORRECTIONASSEMBLER_H
#define BULKEULEREULERPRESSURECORRECTIONASSEMBLER_H

#include "eulerEulerModel.h"
#include "phiAssembler.h"

namespace accel
{

class bulkEulerEulerPressureCorrectionAssembler : public phiAssembler<1>
{
public:
    using Base = phiAssembler<1>;
    using Context = typename Base::Context;
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;

    explicit bulkEulerEulerPressureCorrectionAssembler(
        eulerEulerModel* model);

    // The pressure correction is set up with transportMode `null`, so it has
    // no Gamma to evaluate: skip the base preAssemble, which would otherwise
    // call computeGamma_ and fail on an unset Gamma type. Mirrors
    // pressureCorrectionAssembler in the single-phase solver.
    void preAssemble(const domain* domain, Context* ctx) override
    {
    }

    void postAssemble(const domain* domain, Context* ctx) override;
    void adjustMatrixForPressureReference(const domain* domain, Context* ctx);

protected:
    void assembleNodeTermsFused_(const domain* domain,
                                 Context* ctx) override;
    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;
    void assembleElemTermsInterfaces_(const domain* domain,
                                      Context* ctx) override;
    void assembleElemTermsInterfaceSide_(
        const domain*,
        const interfaceSideInfo*,
        Context*) override
    {
    }
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

private:
    void assembleSideFlux_(const domain* domain,
                           const stk::mesh::PartVector& parts,
                           Context* ctx);
    void assembleSpecifiedPressureCorrection_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx);

    eulerEulerModel* model_;
};

} // namespace accel

#endif // BULKEULEREULERPRESSURECORRECTIONASSEMBLER_H
