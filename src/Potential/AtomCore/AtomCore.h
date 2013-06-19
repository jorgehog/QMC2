/* 
 * File:   AtomCore.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:42 PM
 */

#ifndef ATOMCORE_H
#define	ATOMCORE_H

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

#endif	/* ATOMCORE_H */