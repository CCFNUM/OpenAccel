// File       : postProcessProbe.cpp
// Created    : Tue Aug 05 2025 19:49:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "messager.h"
#include "postProcess.h"

namespace accel
{

probeObject::probeObject(postProcess* postProcessPtr,
                         std::string name,
                         postProcessType type,
                         std::vector<std::string> location,
                         label frequency,
                         bool writeToFile,
                         std::array<scalar, SPATIAL_DIM> probeLocation,
                         std::string field)
    : postProcessObject(postProcessPtr,
                        name,
                        type,
                        location,
                        frequency,
                        writeToFile),
      probeLocation_(probeLocation), field_(field)
{
    errorMsg("Probe is not implemented yet");

    if (messager::master())
    {
        const stk::mesh::MetaData& metaData =
            postProcessPtr_->meshRef().metaDataRef();
        const stk::mesh::BulkData& bulkData =
            postProcessPtr_->meshRef().bulkDataRef();

        // Check if a node field
        const STKScalarField* STKFieldPtr =
            metaData.get_field<scalar>(stk::topology::NODE_RANK, field_);

        if (STKFieldPtr)
        {
            // find encapsulating element
        }
        else
        {
            std::string msg = field_ + " not valid for post-process. "
                                       "Available fields:\n";
            const auto& fields = metaData.get_fields();
            for (auto field : fields)
            {
                if (field->entity_rank() == stk::topology::NODE_RANK)
                {
                    msg += field->name() + "\n";
                }
            }

            errorMsg(msg);
        }

        if (writeToFile_)
        {
            if (fs::is_regular_file(postProcessPtr->directory() / name_))
            {
                fs::remove_all(postProcessPtr->directory() / name_);
            }
            std::string fileName(postProcessPtr->directory() / name_);
            std::ofstream file(fileName);

            file << postProcessPtr_->instanceHeader() + "\t";

            if (STKFieldPtr->max_size() == 1)
            {
                file << "value";
            }
            else
            {
                for (label i = 0; i < STKFieldPtr->max_size(); i++)
                {
                    file << "value_" << std::to_string(i) << "\t";
                }
            }

            file << "\n";
            file.close();
        }
    }
}

void probeObject::update()
{
}

} // namespace accel
