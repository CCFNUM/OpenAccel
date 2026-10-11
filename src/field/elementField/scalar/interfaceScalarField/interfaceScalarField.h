// File       : interfaceScalarField.h
// Created    : Wed Oct 07 2026 13:34:00 (+0300)
// Author     : Mohammad Nasreddine
// Description: Element scalar field that lives on interface/boundary sides
//              only (no interior allocation); holds the Robin-Neumann FSI
//              wall fields: normal acceleration, stored flux Q (Tukovic
//              Eq. 30) and frozen explicit Robin term E
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef INTERFACESCALARFIELD_H
#define INTERFACESCALARFIELD_H

#include "elementScalarField.h"

namespace accel
{

class interfaceScalarField : public elementScalarField
{
public:
    // Constructors

    interfaceScalarField(realm* realmPtr,
                         const std::string name,
                         unsigned numberOfStates);

    // Methods

    void putFieldOnRegisteredParts_() override;
};

} // namespace accel

#endif // INTERFACESCALARFIELD_H
