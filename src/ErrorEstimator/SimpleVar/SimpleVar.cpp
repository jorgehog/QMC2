/* 
 * File:   SimpleVar.cpp
 * Author: jorgmeister
 * 
 * Created on October 30, 2012, 6:13 PM
 */

#include "../../QMCheaders.h"

SimpleVar::SimpleVar(int n_c, ParParams & pp) 
: ErrorEstimator(n_c, "", "", pp.parallel, pp.node, pp.n_nodes, false){
    data_to_file = false;
}



