/* 
 * File:   Brute_Force.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */


#include "Brute_Force.h"

#include "../../structs.h"

#include "../../Diffusion/Simple/Simple.h"


Brute_Force::Brute_Force(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {
    diffusion = new Simple(n_p, dim, 0, gP.random_seed);

}
