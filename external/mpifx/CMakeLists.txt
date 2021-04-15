# This CMakeLists.txt should be invoked via fnet_config_hybrid_dependency()!

add_subdirectory(${MPIFX_SOURCE_DIR} ${MPIFX_BINARY_DIR})

add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)
