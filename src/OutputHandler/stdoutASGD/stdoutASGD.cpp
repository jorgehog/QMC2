/* 
 * File:   stdoutASGD.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:05 PM
 */

#include "../../QMCheaders.h"

stdoutASGD::stdoutASGD(std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int num_procs)
: OutputHandler(filename, path, parallel, my_rank, num_procs) {
    this->is_ASGD = true;

    aGrad = 0;
    bGrad = 0;
    sumE = 0;
}

void stdoutASGD::dump() {
    aGrad += asgd->gradient_tot[0];
    bGrad += asgd->gradient_tot[1];
    sumE += asgd->E;

    file << aGrad / asgd->sample << "\t" << bGrad / asgd->sample << "\t";
    file << fabs(asgd->step) << "\t";
    
    for (int i = 0; i < asgd->Nspatial_params; i++) {
        file << asgd->vmc->get_orbitals_ptr()->get_parameter(i) << "\t";
    }
    for (int i = 0; i < asgd->Njastrow_params; i++){
        file << asgd->vmc->get_jastrow_ptr()->get_parameter(i) << "\t";
    }

    file << asgd->E << "\t";
    file << sumE / asgd->sample << std::endl;
}

