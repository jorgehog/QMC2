/* 
 * File:   Minimizer.cpp
 * Author: jorgehog
 * 
 * Created on 23. august 2012, 16:52
 */

#include "../QMCheaders.h"

Minimizer::Minimizer(VMC* vmc, const arma::rowvec & alpha, const arma::rowvec & beta) {

    this->vmc = vmc;

    this->Nspatial_params = alpha.n_rows;
    this->Njastrow_params = beta.n_rows * vmc->get_jastrow_ptr()->active;
    this->Nparams = this->Njastrow_params + this->Nspatial_params;

    for (int i = 0; i < Nspatial_params; i++) {
        vmc->get_orbitals_ptr()->set_parameter(alpha(i), i);
    }

    for (int i = 0; i < Njastrow_params; i++) {
        vmc->get_jastrow_ptr()->set_parameter(beta(i), i);
    }

}

void Minimizer::output(std::string message, double number) {
    using namespace std;

    if (number != -1) {
        cout << message << " " << number << endl;
    } else {
        cout << message << endl;
    }


    if (Nspatial_params != 0) std::cout << "\nAlpha:\n";
    for (int alpha = 0; alpha < Nspatial_params; alpha++) {
        cout << "\t" << vmc->get_orbitals_ptr()->get_parameter(alpha) << endl;
    }

    if (Njastrow_params != 0) std::cout << "\nBeta:\n";
    for (int beta = 0; beta < Njastrow_params; beta++) {
        cout << "\t" << vmc->get_jastrow_ptr()->get_parameter(beta) << endl;
    }

    cout << endl;

}