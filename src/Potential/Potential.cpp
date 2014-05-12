#include "Potential.h"

using namespace QMC2;

Potential::Potential(int n_p, int dim) {

    this->n_p = n_p;
    this->dim = dim;

    name = "Potential";
    
}

Potential::Potential(){
    
}
