
#ifndef FREEZEOUTSURFACE_H_
#define FREEZEOUTSURFACE_H_

#include <stdlib.h>
#include <vector>

typedef struct
{
    // standard vh format for iS3D (convert to double at end of simulation)

    std::vector<float> tau;                                	// contravariant freezeout cell position
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> eta;

    std::vector<float> dsigma_tau;                         	// covariant surface normal vector
    std::vector<float> dsigma_x;
    std::vector<float> dsigma_y;
    std::vector<float> dsigma_eta;

    std::vector<float> E;                                  	// energy density [GeV/fm^3]
    std::vector<float> T;                                  	// temperature [GeV]
    std::vector<float> P;                                  	// equilibrium pressure [GeV/fm^3]

    std::vector<float> ux;                                 	// contravariant fluid velocity
    std::vector<float> uy;
    std::vector<float> un;
    														// contravariant shear stress pi^{\mu\nu}
    std::vector<float> pixx;                               	// [GeV/fm^3]
    std::vector<float> pixy;                               	// [GeV/fm^3]
    std::vector<float> pixn;                               	// [GeV/fm^4]
    std::vector<float> piyy;                               	// [GeV/fm^3]
    std::vector<float> piyn;                               	// [GeV/fm^4]
    std::vector<float> pinn;                               	// [GeV/fm^5]    (extraneous so leave it empty I guess)

    std::vector<float> Pi;                                 	// bulk pressure [GeV/fm^3]

} freezeout_surface;

#endif







