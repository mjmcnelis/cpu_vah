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
    Regulations.cpp
    SourceTerms.cpp
    TransportAniso.cpp
    TransportAnisoNonconformal.cpp
    TransportViscous.cpp
    Trento.cpp
    Viscosities.cpp
    ViscousBjorken.cpp
    ViscousGubser.cpp
    )

add_library(VAH_lib SHARED ${SOURCES})
add_executable(HydroWrapper ${SOURCES})

target_link_libraries(VAH_lib ${GSL_LIBRARIES})
target_link_libraries(HydroWrapper ${GSL_LIBRARIES})

install(TARGETS HydroWrapper DESTINATION ${CMAKE_HOME_DIRECTORY})
