/* 
 * File:   Harmonic_osc.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:42 PM
 */

#ifndef HARMONIC_OSC_H
#define	HARMONIC_OSC_H

/*! \brief Implementation of the Harmonic Oscillator potential. 0.5*w**2*r**2
 */
class Harmonic_osc : public Potential {
protected:
    
    double w; //!< The oscillator frequency.

public:

    Harmonic_osc(GeneralParams &);

    double get_pot_E(const Walker* walker) const;

};

#endif	/* HARMONIC_OSC_H */

