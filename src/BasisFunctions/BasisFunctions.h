/* 
 * File:   function.h
 * Author: jorgehog
 *
 * Created on 27. juni 2012, 13:17
 */

#ifndef BASISFUNCTIONS_H
#define	BASISFUNCTIONS_H

class BasisFunctions {
public:
    BasisFunctions();
    
    virtual double eval(const Walker* walker, int i) const = 0;
};

#endif	/* BASISFUNCTIONS_H */

