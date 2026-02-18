target_include_directories(main_obj SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/exprtk
    ${CMAKE_CURRENT_SOURCE_DIR}/nanoflann/include
    ${CMAKE_CURRENT_SOURCE_DIR}/gnuplot-cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/include
)

# Nalu master_element sources
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Hex27CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Hex8CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Hex8FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/MasterElement.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/MasterElementFactory.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Pyr5CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Quad42DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Quad43DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Quad92DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tet10FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tet4CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tri32DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tri6FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Wed6CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/MasterElementWork.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateBentHexVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateHexVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateOctohedronVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculatePolyhedralArea.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculatePolyhedralVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateQuadArea.F
)

# Nalu tabular_props sources
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/BSpline.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/Converter.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/Functions.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/H5IO.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/HDF5FilePtr.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/HDF5Table.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/LU.C
)

# Nalu others
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/NaluEnv.C
)
