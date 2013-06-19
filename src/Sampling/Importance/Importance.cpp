/* 
 * File:   Importance.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "Importance.h"

#include "../../structs.h"

#include "../../Diffusion/Fokker_Planck/Fokker_Planck.h"
#include "../../Walker/Walker.h"


Importance::Importance(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {

    diffusion = new Fokker_Planck(n_p, dim, 0, gP.random_seed);

}

void Importance::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
    using namespace arma;

    walker_pre->spatial_grad(span(start, end), span()) = walker_post->spatial_grad(span(start, end), span());
    walker_pre->jast_grad = walker_post->jast_grad;

    walker_pre->qforce = walker_post->qforce;

}

void Importance::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    using namespace arma;

    walker_post->spatial_grad(span(start, end), span()) = walker_pre->spatial_grad(span(start, end), span());
    walker_post->jast_grad = walker_pre->jast_grad;
}

void Importance::copy_walker(const Walker* parent, Walker* child) const {
        child->jast_grad = parent->jast_grad;
        child->spatial_grad = parent->spatial_grad;
        child->qforce = parent->qforce;
    }