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

    return ran3(&random_seed);
}

double Diffusion::get_new_pos(const Walker* walker, int i, int j) {
    return gaussian_deviate(&random_seed) * std;
}

Simple::Simple(int n_p, int dim, double timestep, long random_seed, double D)
: Diffusion(n_p, dim, timestep, random_seed, D) {

}

double Simple::get_new_pos(const Walker* walker, int i, int j) {
    return Diffusion::get_new_pos(walker, i, j);
}

double Simple::get_g_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {
    return 1;
}

double Simple::get_GBfunc(Walker* walker_pre, Walker* walker_post, double E_T) const {
    double Vx = qmc->get_system_ptr()->get_potential_energy(walker_pre);
    double Vy = qmc->get_system_ptr()->get_potential_energy(walker_post);

    return exp(-(0.5 * (Vy + Vx) - E_T) * timestep);


}

Fokker_Planck::Fokker_Planck(int n_p, int dim, double timestep, long random_seed, double D)
: Diffusion(n_p, dim, timestep, random_seed, D) {

}

double Fokker_Planck::get_new_pos(const Walker* walker, int i, int j) {
    return D * timestep * walker->qforce(i, j) + Diffusion::get_new_pos(walker, i, j);
}

double Fokker_Planck::get_g_ratio(const Walker* walker_post, const Walker* walker_pre, int particle) const {

    double g_ratio = 0;
    for (int j = 0; j < dim; j++) {
        g_ratio += 0.5 * (walker_pre->qforce(particle, j) + walker_post->qforce(particle, j))*
                (D * timestep * 0.5 * (walker_pre->qforce(particle, j) - walker_post->qforce(particle, j))
                - walker_post->r(particle, j) + walker_pre->r(particle, j));
    }

    return exp(g_ratio);
}

double Fokker_Planck::get_GBfunc(Walker* walker_pre, Walker* walker_post, double E_T) const {
    
    double Ex = qmc->calculate_local_energy(walker_pre);
    double Ey = qmc->calculate_local_energy(walker_post);

    return exp(-(0.5*(Ey + Ex) - E_T)*timestep*qmc->get_accepted_ratio());
}
