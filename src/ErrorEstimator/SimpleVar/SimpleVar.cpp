/* 
 * File:   SimpleVar.cpp
 * Author: jorgmeister
 * 
 * Created on October 30, 2012, 6:13 PM
 */

#include "../../QMCheaders.h"

SimpleVar::SimpleVar(int n_c, ParParams & pp)
: ErrorEstimator(n_c, "", "", pp.parallel, pp.node, pp.n_nodes, false) {
    data_to_file = false;
}

SimpleVar::SimpleVar(int n_c)
: ErrorEstimator(n_c, "", "", false, 0, 1, false) {
    data_to_file = false;
}

double SimpleVar::estimate_error() {
    if (parallel){
        double var = arma::var(data);
        double mean = arma::mean(data);
        
        return combine_variance(var, mean);
    } else {
        return arma::var(data);
    }
}


