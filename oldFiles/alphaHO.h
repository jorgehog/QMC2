/* 
 * File:   alphaHO.h
 * Author: jorgehog
 *
 * Created on 27. juni 2012, 13:17
 */

#ifndef ALPHAHO_H
#define ALPHAHO_H

class alphaHO : public BasisFunctions {
protected:

    double* k;
    double* k2;

public:

    alphaHO(double* k, double* k2);

};

//START OF GENERATED FUNCTIONS

class alphaHO_1 : public alphaHO {
public:

    alphaHO_1(double* k, double* k2);
    virtual double eval(const Walker* walker, int i) const;

protected:

    double* k;
    double* k2;

};

//END OF GENERATED FUNCTIONS

#endif	/* ALPHAHO_H */