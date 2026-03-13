
add_installed_library(cea_bindc SHARED bindc.F90)
#include(../cmake/CompilerFlags.cmake)
project_enable_fortran_std(cea_bindc)
target_link_libraries(cea_bindc PRIVATE cea::core fbasics::core)
target_include_directories(cea_bindc PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(FILES cea.h cea_enum.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cea/bindc)
include(CheckIncludeFile)

# Optionally enable support for the Fortran C interop header
option(USE_ISO_FORTRAN_BINDING "Check for and enable ISO_Fortran_binding.h support" ON)
if(USE_ISO_FORTRAN_BINDING)
    check_include_file("ISO_Fortran_binding.h" HAVE_ISO_FORTRAN_BINDING_H)
    if(HAVE_ISO_FORTRAN_BINDING_H)
        add_compile_definitions(ISO_Fortran_binding_h)
        message(STATUS "Found ISO_Fortran_binding.h; defining ISO_Fortran_binding_h")
    else()
        message(STATUS "ISO_Fortran_binding.h not found; not defining ISO_Fortran_binding_h")
    endif()
endif()
if(WIN32)
    # on windows, several dll either need to be installed on the target
    # machine, or they need to be statically linked. linking seemed easier than
    # finding the path to the various DLLs and installing them.
    # If the intel compilers are used, static link the intel runtime
    # The appropriate intel flags are different on linux and windows, because

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        target_link_options(cea_bindc
                            PUBLIC "-libs:static")
        set_target_properties(cea_bindc PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS ON
            PREFIX ""
            MSVC_RUNTIME_LIBRARY "MultiThreaded"
        )
    else()
        target_link_options(cea_bindc
                            PUBLIC "-static-libgfortran"
                            PUBLIC "-static-libgcc"
                            PUBLIC "-static-libquadmath")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_link_options(cea_bindc
                        PUBLIC "-static-intel")
endif()

if (CEA_BUILD_TESTING)

    add_executable(cea_bindc_rp1311_ex1 samples/rp1311_example1.c)
    target_link_libraries(cea_bindc_rp1311_ex1 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex1
        COMMAND cea_bindc_rp1311_ex1
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

        if(HAVE_ISO_FORTRAN_BINDING_H)
            add_executable(cea_bindc_rp1311_ex1_test_structure samples/rp1311_example1_test_structure.c)
            target_link_libraries(cea_bindc_rp1311_ex1_test_structure PRIVATE cea::bindc)
            add_test(
                    NAME cea_bindc_rp1311_ex1_test_structure
                    COMMAND cea_bindc_rp1311_ex1_test_structure
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            )
        else()
            message(STATUS "Skipping cea_bindc_rp1311_ex1_test_structure: ISO_Fortran_binding.h not found")
        endif()

    add_executable(cea_bindc_rp1311_ex2 samples/rp1311_example2.c)
    target_link_libraries(cea_bindc_rp1311_ex2 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex2
        COMMAND cea_bindc_rp1311_ex2
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(cea_bindc_rp1311_ex3 samples/rp1311_example3.c)
    target_link_libraries(cea_bindc_rp1311_ex3 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex3
        COMMAND cea_bindc_rp1311_ex3
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(cea_bindc_rp1311_ex5 samples/rp1311_example5.c)
    target_link_libraries(cea_bindc_rp1311_ex5 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex5
        COMMAND cea_bindc_rp1311_ex5
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(cea_bindc_rp1311_ex6 samples/rp1311_example6.c)
    target_link_libraries(cea_bindc_rp1311_ex6 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex6
        COMMAND cea_bindc_rp1311_ex6
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(cea_bindc_rp1311_ex7 samples/rp1311_example7.c)
    target_link_libraries(cea_bindc_rp1311_ex7 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex7
        COMMAND cea_bindc_rp1311_ex7
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(cea_bindc_rp1311_ex8 samples/rp1311_example8.c)
    target_link_libraries(cea_bindc_rp1311_ex8 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex8
        COMMAND cea_bindc_rp1311_ex8
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    add_executable(cea_bindc_rp1311_ex9 samples/rp1311_example9.c)
    target_link_libraries(cea_bindc_rp1311_ex9 PRIVATE cea::bindc)
    add_test(
        NAME cea_bindc_rp1311_ex9
        COMMAND cea_bindc_rp1311_ex9
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

endif()
