#pragma once


#include "../Potential.h"


namespace QMC2
{


struct GeneralParams;

class DiAtomCore : public Potential {
public:
    DiAtomCore(GeneralParams & gp);
    
    double get_pot_E(const Walker* walker) const;
    
private:
    
    double *R; //<! Distance between cores.
    int Z;
    
};

}
