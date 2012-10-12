/* 
 * File:   Jastrow.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 16:52
 */

#include "../QMCheaders.h"


Jastrow::Jastrow(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = n_p / 2;
    this->dim = dim;
}

Jastrow::Jastrow() {

}
