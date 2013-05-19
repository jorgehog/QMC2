/* 
 * File:   DoubleWell.cpp
 * Author: jorgmeister
 * 
 * Created on May 10, 2013, 1:25 PM
 */

#include "../../QMCheaders.h"

DoubleWell::DoubleWell(GeneralParams& gp) :
Potential(gp.n_p, gp.dim){
    
    this->w = gp.systemConstant; 
    this->R = &(gp.R);
    
    name = "DoubleWell";
    
}


double DoubleWell::get_pot_E(const Walker* walker) const {
    
    double e_pot = 0;
    
    double quarterR2 = 0.25*(*R)*(*R);
    
    for (int i = 0; i < n_p; i++) {
        
        e_pot += walker->get_r_i2(i)+ quarterR2 - (*R)*fabs(walker->r(i, 0));

    }

    
    return 0.5*w*w*e_pot;
    
}               