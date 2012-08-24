/* 
 * File:   MINIMIZER.h
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#ifndef MINIMIZER_H
#define	MINIMIZER_H

class Minimizer {
protected:
    
    VMC* vmc;
    
    int Nspatial_params;
    int Njastrow_params;
    int Nparams;
    
    //virtual void get_variational_gradients(); is this a general req?
    
public:
    
    Minimizer(VMC* vmc, int NSP, int NJP);
    
    Orbitals* get_orbitals(){
        return vmc->get_orbitals_ptr();
    }
    
    Jastrow* get_jastrow() {
        return vmc->get_jastrow_ptr();
    }
    
    virtual VMC* minimize() = 0;
    
    void output();
    
};


#endif	/* MINIMIZER_H */