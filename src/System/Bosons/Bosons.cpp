/* 
 * File:   Bosons.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:44 PM
 */

#include "../../QMCheaders.h"
#include "Bosons.h"

Bosons::Bosons(GeneralParams & gP, Orbitals* orbital)
: System(gP.n_p, gP.dim, orbital) {

    overlap = false;

}

void Bosons::get_spatial_grad_full(Walker* walker) const {

    for (int i = 0; i < n_p; i++) {
        for (int k = 0; k < dim; k++) {
            walker->spatial_grad(i, k) = walker->dell_phi(i)(k) / walker->phi(i);
        }
    }

}

void Bosons::get_spatial_grad(Walker* walker, int particle) const {

    for (int k = 0; k < dim; k++) {
        walker->spatial_grad(particle, k) = walker->dell_phi(particle)(k) / walker->phi(particle);
    }

}

double Bosons::get_spatial_wf(const Walker* walker) {

    double wf = 1;

    //Using the phi matrix as a vector in the case of bosons.
    //Assuming all particles to occupy the same single particle state (neglecting permutations).
    for (int i = 0; i < n_p; i++) {
        wf *= walker->phi(i);
    }

    return wf;

}

double Bosons::get_spatial_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) {
    double s_ratio;

    //Assumes every single particle state to be identical.
    s_ratio = walker_post->phi(particle) / walker_pre->phi(particle);

    a = false; //FIX

    return s_ratio;
}

double Bosons::get_spatial_lapl_sum(Walker* walker) const {

    double sum = 0;
    for (int i = 0; i < n_p; i++) {
        sum += orbital->lapl_phi(walker, i, i) / walker->phi(i);
    }

    return sum;
}
