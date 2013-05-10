/* 
 * File:   DiAtomCore.h
 * Author: jorgmeister
 *
 * Created on May 10, 2013, 1:25 PM
 */

#ifndef DIATOMCORE_H
#define	DIATOMCORE_H

class DiAtomCore : public Potential {
public:
    DiAtomCore(GeneralParams & gp, double* R);
    
    double get_pot_E(const Walker* walker) const;
    
private:
    
    double *R; //<! Distance between cores.
    int Z;
    
};

#endif	/* DIATOMCORE_H */

