/* 
 * File:   DoubleWell.h
 * Author: jorgmeister
 *
 * Created on May 10, 2013, 1:25 PM
 */

#ifndef DOUBLEWELL_H
#define	DOUBLEWELL_H

#include "../Potential.h"
struct GeneralParams;

class DoubleWell : public Potential {
public:
    DoubleWell(GeneralParams & gp);
    
    double get_pot_E(const Walker* walker) const;
    
private:
 
    double *R; //<! Distance between wells.
    double w; //< Oscillator frequency
    
};

#endif	/* DOUBLEWELL_H */