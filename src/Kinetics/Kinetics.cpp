/* 
 * File:   Kinetics.cpp
 * Author: jorgehog
 * 
 * Created on 13. april 2012, 17:45
 */

#include "../QMCheaders.h"


Kinetics::Kinetics() {

}

Kinetics::Kinetics(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = n_p / 2;
    this->dim = dim;
}

