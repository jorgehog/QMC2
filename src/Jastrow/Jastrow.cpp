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

void Jastrow::get_dJ_matrix(Walker* walker) const {
    for (int i = 0; i < n_p; i++) {
        this->get_dJ_matrix(walker, i);
    }
}

double Jastrow::get_variational_derivative(const Walker* walker, int n){
    double b = get_parameter(n);
    double h = 0.0001;
    
    set_parameter(b + h, n);
    
    double phip = get_val(walker);
    set_parameter(b - h, n);
    
    double phim = get_val(walker);
    set_parameter(b, n);
    
    return (phip - phim) / (2 * h * get_val(walker));
}