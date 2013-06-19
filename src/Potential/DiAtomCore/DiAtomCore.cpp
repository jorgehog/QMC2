/* 
 * File:   DiAtomCore.cpp
 * Author: jorgmeister
 * 
 * Created on May 10, 2013, 1:25 PM
 */

#include "DiAtomCore.h"

#include "../../misc.h"
#include "../../Walker/Walker.h"


DiAtomCore::DiAtomCore(GeneralParams& gp) :
Potential(gp.n_p, gp.dim){
    
    this->Z = gp.n_p/2; 
    this->R = &(gp.R);
    
    name = "DiAtomCore";
    
}


double DiAtomCore::get_pot_E(const Walker* walker) const {
    
    double e_pot = 0;
    double com_corr, shared;
    
    double quarterR2 = 0.25*(*R)*(*R);
    for (int i = 0; i < n_p; i++) {
        
        shared = walker->get_r_i2(i)+ quarterR2;
        com_corr = (*R)*walker->r(i, 0);
        
        e_pot -= Z*(1./sqrt(shared + com_corr) + 1./sqrt(shared - com_corr));
//        e_pot -= 2*Z*shared/sqrt(shared*shared - com_corr*com_corr);
    }

    e_pot += Z*Z/(*R);
    
    return e_pot;
    
}
