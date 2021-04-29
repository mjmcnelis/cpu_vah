set (SOURCES
    AdaptiveTimeStep.cpp
    AnisoBjorken.cpp
    AnisoGubser.cpp
    AnisoVariables.cpp
    DynamicalVariables.cpp
    EquationOfState.cpp
    FluxTerms.cpp
    FreezeoutFinder.cpp
    GhostCells.cpp
    Hydrodynamics.cpp
    HydroWrapper.cpp
    IdealGubser.cpp
    InferredVariables.cpp
    InitialConditions.cpp
    KurganovTadmor.cpp
    Main.cpp
    Memory.cpp
    NeighborCells.cpp
    Output.cpp
    Parameters.cpp
    Print.cpp
    Projections.cpp
    Regulation.cpp
    SourceTerms.cpp
    TransportAniso.cpp
    TransportAnisoNonconformal.cpp
    TransportViscous.cpp
    Trento.cpp
    Viscosities.cpp
    ViscousBjorken.cpp
    ViscousGubser.cpp
    )

add_library(cpu_vah_lib SHARED ${SOURCES})
add_executable(cpu_vah ${SOURCES})

target_link_libraries(cpu_vah_lib ${GSL_LIBRARIES})
target_link_libraries(cpu_vah ${GSL_LIBRARIES})

install(TARGETS cpu_vah DESTINATION ${CMAKE_HOME_DIRECTORY})
