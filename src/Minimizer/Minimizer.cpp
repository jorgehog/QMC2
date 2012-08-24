/* 
 * File:   Minimizer.cpp
 * Author: jorgehog
 * 
 * Created on 23. august 2012, 16:52
 */

#include "../QMCheaders.h"

Minimizer::Minimizer(VMC* vmc, int NSP, int NJP) {
    this->vmc = vmc;
    this->Njastrow_params = NJP;
    this->Nspatial_params = NSP;
    this->Nparams = NJP + NSP;
}

void Minimizer::output() {
    std::cout << "Finished minimizing. Final parameters:";
    
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