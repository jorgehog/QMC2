/* 
 * File:   SimpleVar.cpp
 * Author: jorgmeister
 * 
 * Created on October 30, 2012, 6:13 PM
 */

#include "SimpleVar.h"

#include "../../structs.h"

SimpleVar::SimpleVar(ParParams & pp)
: ErrorEstimator(0, "", "", pp.parallel, pp.node, pp.n_nodes, false) {
    data_to_file = false;
    f = 0;
    f2 = 0;
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

void SimpleVar::normalize(){
    
}