#include "Potential.h"

using namespace QMC2;

Potential::Potential(int n_p, int dim, std::string name) :
    name(name)
{

    this->n_p = n_p;
    this->dim = dim;

    pot_sampler = new Sampler(name);
    
}

Potential::Potential(){
    
}

std::string Potential::name_suffix;
