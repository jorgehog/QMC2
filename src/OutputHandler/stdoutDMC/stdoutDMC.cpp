/* 
 * File:   stdoutDMC.cpp
 * Author: jorgmeister
 * 
 * Created on September 4, 2012, 7:22 PM
 */

#include "../../QMCheaders.h"

stdoutDMC::stdoutDMC(std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int n_nodes)
: OutputHandler(filename, path, parallel, my_rank, n_nodes) {

    this->is_dmc = true;
    sumE = 0;
    sumN = 0;
    n = 0;

}

void stdoutDMC::dump() {
    
    sumE += dmc->dmc_E / dmc->cycle;
    sumN += dmc->n_w;
    n++;

    file << dmc->dmc_E / dmc->cycle << "\t";
    file << sumE / n << "\t";
    file << dmc->n_w << "\t";
    file << sumN / n << "\t";
    file << dmc->E_T << std::endl;
    

}