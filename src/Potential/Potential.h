/* 
 * File:   Potential.h
 * Author: jorgehog
 *
 * Created on 13. april 2012, 17:21
 */

#ifndef POTENTIAL_H
#define	POTENTIAL_H

class Potential {
protected:
    int n_p;
    int dim;

public:
    Potential(int n_p, int dim);
    Potential();

    virtual double get_pot_E(const Walker* walker) const = 0;

};

#endif	/* POTENTIAL_H */

