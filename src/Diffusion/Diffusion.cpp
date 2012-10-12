/* 
 * File:   Diffusion.cpp
 * Author: jorgehog
 * 
 * Created on 16. april 2012, 14:03
 */

#include "../QMCheaders.h"


Diffusion::Diffusion(int n_p, int dim, double timestep, long random_seed, double D) {
    this->n_p = n_p;
    this->dim = dim;
    this->timestep = timestep;
    this->random_seed = random_seed;
    this->D = D;
    this->std = sqrt(2 * D * timestep);
}

double Diffusion::call_RNG() {
    return ran2(&random_seed);
}

double Diffusion::get_new_pos(const Walker* walker, int i, int j) {
    return gaussian_deviate(&random_seed) * std;
}

double Diffusion::get_GBfunc(double E_x, double E_y, double E_T) const {
    return exp(-(0.5 * (E_x + E_y) - E_T) * timestep);
}



