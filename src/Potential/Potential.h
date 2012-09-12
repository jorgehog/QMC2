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

    virtual double get_pot_E(const Walker* walker) const = 0;

};

class Harmonic_osc : public Potential {
protected:
    double w;

public:

    Harmonic_osc(int n_p, int dim, double W);

    virtual double get_pot_E(const Walker* walker) const;

};

class Coulomb : public Potential {
public:

    Coulomb(int n_p, int dim);

    virtual double get_pot_E(const Walker* walker) const;

};


#endif	/* POTENTIAL_H */

