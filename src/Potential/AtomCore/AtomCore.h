/* 
 * File:   AtomCore.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:42 PM
 */

#ifndef ATOMCORE_H
#define	ATOMCORE_H

class AtomCore : public Potential {
protected:
    int Z;

public:

    AtomCore(GeneralParams &);

    virtual double get_pot_E(const Walker* walker) const;

};

#endif	/* ATOMCORE_H */