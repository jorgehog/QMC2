/* 
 * File:   Coulomb.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:41 PM
 */

#ifndef COULOMB_H
#define	COULOMB_H

class Coulomb : public Potential {
public:

    Coulomb(GeneralParams &);

    virtual double get_pot_E(const Walker* walker) const;

};

#endif	/* COULOMB_H */

