#pragma once


#include "../Potential.h"

struct GeneralParams;

/*! \brief Implementation of the Atom Core potential. -Z/r
 */
class AtomCore : public Potential {
protected:
    int Z; //!< The core charge.
    
public:

    AtomCore(GeneralParams &);

    double get_pot_E(const Walker* walker) const;

};
