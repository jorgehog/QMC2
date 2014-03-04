#pragma once


#include "../Potential.h"


namespace QMC2
{


struct GeneralParams;

class DoubleWell : public Potential {
public:
    DoubleWell(GeneralParams & gp, double R);
    
    double get_pot_E(const Walker* walker) const;
    
private:
 
    double *R; //<! Distance between wells.
    double w; //< Oscillator frequency
    
};

}
