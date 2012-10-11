/* 
 * File:   alphaHO.cpp
 * Author: jorgehog
 * 
 * Created on 27. juni 2012, 13:17
 */

#include "../../QMCheaders.h"

alphaHO::alphaHO(double* k, double* k2) {
    this->k = k;
    this->k2 = k2;
}


//START OF CONSTRUCTORS

alphaHO_1::alphaHO_1(double* k, double* k2)
: alphaHO(k, k2) {

}

//END OF CONSTRUCTORS
//START OF EVAL FUNCTIONS

double alphaHO_1::eval(const Walker* walker, int i) const {
    return 0;
}

//END OF EVAL FUNCTIONS