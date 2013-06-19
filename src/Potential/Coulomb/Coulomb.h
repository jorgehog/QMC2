/* 
 * File:   Coulomb.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:41 PM
 */

#ifndef COULOMB_H
#define	COULOMB_H

#include "../Potential.h"
struct GeneralParams;

/*! \brief Implementation of the Coulomb potential. 1/r_{ij}
 */
class Coulomb : public Potential {
public:

    Coulomb(GeneralParams &);

    double get_pot_E(const Walker* walker) const;

};

#endif	/* COULOMB_H */

