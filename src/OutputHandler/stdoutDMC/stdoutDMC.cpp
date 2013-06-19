/* 
 * File:   stdoutDMC.cpp
 * Author: jorgmeister
 * 
 * Created on September 4, 2012, 7:22 PM
 */

#include "stdoutDMC.h"

#include "../../QMC/DMC/DMC.h"


stdoutDMC::stdoutDMC(std::string path, std::string filename)
: OutputHandler(filename, path, false, 0, 1) {

    this->is_dmc = true;
    sumE = 0;
    sumN = 0;
    n = 0;
    
    init_file();

}

void stdoutDMC::dump() {
    
    sumE += dmc->dmc_E;
    sumN += dmc->n_w_tot;
    n++;

    file << dmc->dmc_E << "\t";
    file << sumE / n << "\t";
    file << dmc->n_w_tot << "\t";
    file << sumN / n << "\t";
    file << dmc->E_T << std::endl;

}