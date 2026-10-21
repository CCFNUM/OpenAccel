// File       : fsiConvergence.h
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: FSI convergence checks
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef FSICONVERGENCE_H
#define FSICONVERGENCE_H

#include "types.h"

namespace accel
{

class simulation;

class fsiConvergence
{
public:
    explicit fsiConvergence(simulation& sim);

    bool enabled() const;

    void resetForTimeStep();

    void update();

    // FSI interface residual of the solid predictor (iteration 0): enters the
    // timestep maxima of the dual norm (solids4foam updateResidual)
    void updatePredictor();

    // false: convergence falls back to the equation residuals
    bool hasCriteria() const;

    bool isConverged() const;

private:
    simulation& sim_;

    std::map<label, scalar> fsiInterfaceResidualNormMax_;
    // running max of ||D|| (second norm of the dual residual)
    std::map<label, scalar> fsiInterfaceMaxTotalDisplNorm_;
    std::map<label, scalar> fsiInterfaceResidualNorms_;
    scalar fsiInterfaceResidualNorm_ = 0.0;

    std::map<label, std::vector<scalar>> fsiInterfaceTractionPrev_;
    std::map<label, scalar> fsiInterfaceTractionResidualNormMax_;
    std::map<label, scalar> fsiForceResidualNorms_;
    scalar fsiForceResidualNorm_ = 0.0;

    std::map<std::string, std::map<label, std::shared_ptr<std::ofstream>>>
        residualStreams_;

    // fsi_diagnostics.out writer
    std::shared_ptr<std::ofstream> diagnosticsStream_;
    bool diagnosticsEnabled_ = false;

    // any fluid-solid interface uses Robin-Neumann (cached)
    bool isRobinNeumann_ = false;

    // RN interface with invalid parameters: DN fallback (cached)
    bool rnFallback_ = false;

    // p' change at the interface, normalised at k=2 (fsi_diagnostics.out)
    std::map<label, std::vector<scalar>> fsiInterfacePCorrPrev_;
    scalar interfacePCorrChangeNorm_ = 0.0;
    scalar interfacePCorrChangeNormRef_ = 1.0; // set on k=2 of each timestep

    // RMS of rho_s*h_s*a_n over the RN BIPs (fsi_diagnostics.out)
    scalar robinAccelTermRms_ = 0.0;

    // RMS of Q_vol/|A| over the RN BIPs (fsi_diagnostics.out)
    scalar robinPCorrFluxRms_ = 0.0;

    // any interface requests write_force_balance
    bool writeForceBalance_ = false;

    // fsi_force_balance.out writer
    std::shared_ptr<std::ofstream> forceBalanceStream_;

    void updateFsiInterfaceResidual_(bool writeResiduals);

    void updateInterfacePCorrChange_();

    void updateFsiForceResidual_(bool writeResiduals);

    // global fluid/solid interface force balance, one row per outer iteration
    void updateForceBalance_();

    void initializeForceBalanceFile_();

    void writeForceBalanceRow_(const std::vector<scalar>& fluidForce,
                               const std::vector<scalar>& solidForce);

    void initializeResidualFile_(label interfIdx,
                                 const std::string& interfName,
                                 const std::string& criterionName);

    void writeResidualLine_(const std::string& criterionName,
                            label interfIdx,
                            scalar residualNorm);

    void initializeDiagnosticsFile_();

    void writeDiagnosticsRow_(label fsiIter,
                              scalar kinematicRes,
                              scalar forceRes,
                              scalar interfacePCorrChange,
                              scalar robinAccelTermRms,
                              scalar robinPCorrFluxRms);
};

} // namespace accel

#endif // FSICONVERGENCE_H
