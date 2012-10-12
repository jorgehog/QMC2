/* 
 * File:   Simple.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:38 PM
 */

#include "../../QMCheaders.h"

Simple::Simple(int n_p, int dim, double timestep, long random_seed, double D)
: Diffusion(n_p, dim, timestep, random_seed, D) {

}

double Simple::get_new_pos(const Walker* walker, int i, int j) {
    return Diffusion::get_new_pos(walker, i, j);
}

double Simple::get_g_ratio(const Walker* walker_post, const Walker* walker_pre) const {
    return 1;
}



