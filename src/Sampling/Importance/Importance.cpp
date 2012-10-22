/* 
 * File:   Importance.cpp
 * Author: jorgmeister
 * 
 * Created on October 12, 2012, 2:43 PM
 */

#include "../../QMCheaders.h"

Importance::Importance(GeneralParams & gP)
: Sampling(gP.n_p, gP.dim) {
    diffusion = new Fokker_Planck(n_p, dim, 0, gP.random_seed, gP.D);

}

void Importance::update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
    using namespace arma;

    int start = n2 * (particle >= n2);
    int end = start + n2 - 1;

    walker_pre->spatial_grad(span(start, end), span()) = walker_post->spatial_grad(span(start, end), span());
    walker_pre->jast_grad = walker_post->jast_grad;

    walker_pre->qforce = walker_post->qforce;

}

void Importance::reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    using namespace arma;

    int start = n2 * (particle >= n2);
    int end = start + n2 - 1;

    walker_post->spatial_grad(span(start, end), span()) = walker_pre->spatial_grad(span(start, end), span());
    walker_post->jast_grad = walker_pre->jast_grad;
}
