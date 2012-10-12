/* 
 * File:   Harmonic_osc.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:42 PM
 */

#ifndef HARMONIC_OSC_H
#define	HARMONIC_OSC_H

class Harmonic_osc : public Potential {
protected:
    double w;

public:

    Harmonic_osc(GeneralParams &);

    virtual double get_pot_E(const Walker* walker) const;

};

#endif	/* HARMONIC_OSC_H */

