#include "Minimizer.h"

#include "../structs.h"

#include "../Orbitals/Orbitals.h"
#include "../Jastrow/Jastrow.h"

using namespace QMC2;

void Minimizer::initializeParameters()
{
    for (int i = 0; i < Nspatial_params; i++) {
        vmc->get_orbitals_ptr()->set_parameter(alpha0(i), i);
    }

    for (int i = 0; i < Njastrow_params; i++) {
        vmc->get_jastrow_ptr()->set_parameter(beta0(i), i);
    }
}

Minimizer::Minimizer(VMC* vmc, const ParParams & pp, const arma::rowvec & alpha0, const arma::rowvec & beta0) {

    this->vmc = vmc;

    this->Nspatial_params = alpha0.n_elem;
    this->Njastrow_params = beta0.n_elem * vmc->get_jastrow_ptr()->active;
    this->Nparams = this->Njastrow_params + this->Nspatial_params;

    this->alpha0 = alpha0;
    this->beta0 = beta0;

    is_master = pp.is_master;
    n_nodes = pp.n_nodes;

}
