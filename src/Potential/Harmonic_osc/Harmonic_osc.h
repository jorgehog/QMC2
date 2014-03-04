#pragma once


#include "../Potential.h"
struct GeneralParams;

/*! \brief Implementation of the Harmonic Oscillator potential. 0.5*w**2*r**2
 */
class Harmonic_osc : public Potential {
protected:
    
    double w; //!< The oscillator frequency.

public:

    Harmonic_osc(GeneralParams &);

    void set_w(double w)
    {
        this->w = w;
    }

    double get_pot_E(const Walker* walker) const;

};

