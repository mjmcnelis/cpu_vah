
#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_

#include "DynamicalVariables.h"
#include "Parameters.h"
#include "Macros.h"

// Butcher tableau of explicit Runge Kutta scheme with s stages (below r = s - 1)

//  c1  |   1   a11
//  c2  |   a20 a21 a22
//  ..  |   ..  ..  ..  ..
//  cr  |   ar0 ar1 ..  ..  arr 0        
//_________________________
//      |   b0  b1  b2  ..  ..  bs

// for the RK update, we take a weighted combination of the q-arguments of each stage (except the final term is the last intermediate Euler step)
// q_{n+1} = b0.q + b1.q1 + ... + bs.dt.E(t + cr.dt, qr)

// this is equivalent to the usual k_n stage notation

#ifdef SHU_OSHER                // Shu-Osher: third-order, strong stability preserving (SSP), 3 stages
    
    #define STAGES 3

    #define C1 1.
    #define C2 1./2.

    #define A11 1.

    #define A20 3./4.
    #define A21 1./4.
    #define A22 1./4.

    #define B0 1./3.
    #define B1 0.
    #define B2 2./3.
    #define B3 2./3.

#elif defined(SPITERI_RUUTH)    // Spiteri-Ruuth: third-order, SSP, 4 stages
    #define STAGES 4

    #define C1 1./2.
    #define C2 1.
    #define C3 1./2.

    #define A11 1./2.

    #define A20 0.
    #define A21 1.
    #define A22 1./2.

    #define A30 2./3.
    #define A31 0.
    #define A32 1./3.
    #define A33 1./6.

    #define B0 0.
    #define B1 0.
    #define B2 0.
    #define B3 1.
    #define B4 1./2.

#else                           // default is Heun: second-order, SSP, 2 stages
    #define HEUN                // used to skip allocation of uI variable
    #define STAGES 2

    #define C1 1.

    #define A11 1.

    #define B0 1./2.
    #define B1 1./2.
    #define B2 1./2.
#endif

void get_RK_stage_1(const hydro_variables * const __restrict__ q, hydro_variables * const __restrict__ q1, precision a1, lattice_parameters lattice);

void get_RK_stage_2(const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ q1, hydro_variables * const __restrict__ q2, precision a0, precision a1, precision a2, lattice_parameters lattice);

void get_RK_stage_3(const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ q1, const hydro_variables * const __restrict__ q2, hydro_variables * const __restrict__ q3, precision a0, precision a1, precision a2, precision a3, lattice_parameters lattice);

void get_RK_stage_4(const hydro_variables * const __restrict__ q, const hydro_variables * const __restrict__ q1, const hydro_variables * const __restrict__ q2, const hydro_variables * const __restrict__ q3, hydro_variables * const __restrict__ q4, precision a0, precision a1, precision a2, precision a3, precision a4, lattice_parameters lattice);

#endif