// File : initialConditions.h
// Created : Wed Mar 06 2024 16:19:10 (+0100)
// Author : Fabian Wermelinger
// Description: Initial conditions utilities
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H

// code
#include "domain.h"
#include "macros.h"
#include "nodeField.h"
#include "types.h"

namespace accel
{
namespace initialCondition
{

// NOTE: Temporary helper — initial
// conditions should be set in model classes without needing to call
// field.setupInitializationDetails(...) method for fields

template <class TField>
void setupFieldInitializationOverDomainFromConfig(TField& field,
                                                  const label domain_index,
                                                  const YAML::Node ic_config)
{
    // get ref to initial condition dict of field
    initialConditionDictionary& initCond =
        field.initialConditionRef(domain_index);

    initCond.setType(convertInitialConditionOptionFromString(
        ic_config["option"].template as<std::string>()));

    switch (initCond.type())
    {
        case initialConditionOption::value:
            {
                // extract primary field name (i.e. exlude any extensions)
                std::size_t pos = field.name().find('.');

                std::string fieldName;
                if (pos != std::string::npos)
                {
                    fieldName = field.name().substr(0, pos);
                }
                else
                {
                    fieldName = field.name();
                }

                // query
                initCond.query<TField::NComponents>(
                    ic_config, field.name(), fieldName);
            }
            break;

        case initialConditionOption::null:
            break;
    }
}

template <class TField>
void setupFieldInitializationOverDomainFromValues(TField& field,
                                                  const label domain_index,
                                                  std::vector<scalar> value)
{
    initialConditionDictionary& initCond =
        field.initialConditionRef(domain_index);
    initCond.setType(initialConditionOption::value);
    initCond.setConstantValue<TField::NComponents>(field.name(), value);
}

template <class TField>
void setupFieldInitializationOverDomainFromInput(
    TField& field,
    const std::string& field_name,
    const std::shared_ptr<domain>& domain)
{
    const YAML::Node field_conf =
        domain->getYAMLInitialConditions()[field_name];
    if (field_conf)
    {
        initialCondition::setupFieldInitializationOverDomainFromConfig(
            field, domain->index(), field_conf);
    }
    else
    {
        errorMsg("initialCondition::"
                 "setupFieldInitializationOverDomainFromInput: initialization "
                 "config does not "
                 "define `" +
                 field_name + "` mapping for domain `" + domain->name() + "`");
    }
}

template <class TField>
void setupFluidSpecificFieldInitializationOverDomainFromInput(
    TField& field,
    const std::string& field_name,
    const std::string& phase_name,
    const std::shared_ptr<domain>& domain)
{
    const YAML::Node phase_conf =
        domain->getYAMLInitialConditions()["fluid_specific_initialization"]
                                          [phase_name];
    const YAML::Node field_conf = phase_conf[field_name];
    if (field_conf)
    {
        initialCondition::setupFieldInitializationOverDomainFromConfig(
            field, domain->index(), field_conf);
    }
    else
    {
        errorMsg("initialCondition::"
                 "setupFluidSpecificFieldInitializationOverDomainFromInput: "
                 "initialization "
                 "config does not "
                 "define `" +
                 field_name + "` mapping for domain `" + domain->name() + "`");
    }
}

} /* namespace initialCondition */
} /* namespace accel */

#endif // INITIALCONDITIONS_H
