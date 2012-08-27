/* 
 * File:   Minimizer.cpp
 * Author: jorgehog
 * 
 * Created on 23. august 2012, 16:52
 */

#include "../QMCheaders.h"

Minimizer::Minimizer(VMC* vmc, const rowvec & alpha, const rowvec & beta,
        int NSP, int NJP) {
    this->vmc = vmc;
    this->Njastrow_params = NJP;
    this->Nspatial_params = NSP;
    this->Nparams = NJP + NSP;

    for (int i = 0; i < Nspatial_params; i++) {
        vmc->get_orbitals_ptr()->set_parameter(alpha(i), i);
    }

    for (int i = 0; i < Njastrow_params; i++) {
        vmc->get_jastrow_ptr()->set_parameter(beta(i), i);
    }
    
    cout << vmc->get_orbitals_ptr()->get_parameter(0) << endl;

    
}

void Minimizer::output(std::string message, double number) {
    
    if (number != -1) {
        cout << message << " " << number << endl;
    } else {
        cout << message << endl;
    }


    std::cout << "\nAlpha:\n";
    for (int alpha = 0; alpha < Nspatial_params; alpha++) {
        cout << "\t" << vmc->get_orbitals_ptr()->get_parameter(alpha) << std::endl;
    }

    std::cout << "\nBeta:\n";
    for (int beta = 0; beta < Njastrow_params; beta++) {
        cout << "\t" << vmc->get_jastrow_ptr()->get_parameter(beta) << std::endl;
    }

    std::cout << std::endl;

}