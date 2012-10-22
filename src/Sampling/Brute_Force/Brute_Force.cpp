/* 
 * File:   Brute_Force.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"


Brute_Force::Brute_Force(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {
    diffusion = new Simple(n_p, dim, 0, gP.random_seed, gP.D);

}
