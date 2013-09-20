/* 
 * File:   stdoutASGD.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:05 PM
 */

#include "stdoutASGD.h"

#include "../../Minimizer/ASGD/ASGD.h"
#include "../../Orbitals/Orbitals.h"
#include "../../Jastrow/Jastrow.h"

stdoutASGD::stdoutASGD(ASGD* asgd, std::string path)
: OutputHandler("ASGD_out", path, false, 0, 1) {

    this->asgd = asgd;
    grad = arma::zeros(asgd->Nparams);
    
    sumE = 0;

    init_file();
}

void stdoutASGD::dump() {

    file << asgd->E << "\t";

    sumE += asgd->E;
    file << sumE / asgd->sample << "\t";


    file << fabs(asgd->step);

    for (int i = 0; i < asgd->Nspatial_params; i++) {
        file << "\t" << asgd->vmc->get_orbitals_ptr()->get_parameter(i) << "\t";
        grad(i) += asgd->gradient_tot(i);
        file << grad(i) / asgd->sample;
    }
    for (int i = 0; i < asgd->Njastrow_params; i++) {
        file << "\t" << asgd->vmc->get_jastrow_ptr()->get_parameter(i) << "\t";
        grad(i + asgd->Nspatial_params) += asgd->gradient_tot(i + asgd->Nspatial_params);
        file << grad(i + asgd->Nspatial_params) / asgd->sample;
    }
    file << std::endl;


}

void stdoutASGD::reset()
{
    grad = arma::zeros(asgd->Nparams);

    sumE = 0;

    OutputHandler::reset();
}



