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
    f = 0;
    f2 = 0;
}

SimpleVar::SimpleVar(int n_c)
: ErrorEstimator(n_c, "", "", false, 0, 1, false) {
    data_to_file = false;
}

double SimpleVar::estimate_error() {
    
    //sample variance
    double mean = f/i;
    double var = f2/(i-1) - i*mean*mean/(i-1);
    
    return combine_variance(var, mean, i);
     
}


void SimpleVar::update_data(double val) {
    f2 += val*val;
    f += val;
    i++;
}