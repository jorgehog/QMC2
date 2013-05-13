/* 
 * File:   Jastrow.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 16:52
 */

#include "../QMCheaders.h"

Jastrow::Jastrow(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = ceil(n_p / 2.0);
    this->dim = dim;
}

Jastrow::Jastrow() {

}

void Jastrow::get_dJ_matrix(Walker* walker) const {
    for (int i = 0; i < n_p; i++) {
        this->get_dJ_matrix(walker, i);
    }
}

double Jastrow::get_variational_derivative(const Walker* walker, int n) {
    double b = get_parameter(n);
    double h = 0.0001;

    set_parameter(b + h, n);

    double phip = get_val(walker);
    set_parameter(b - h, n);

    double phim = get_val(walker);
    set_parameter(b, n);

    return (phip - phim) / (2 * h * get_val(walker));
}

double Jastrow::get_derivative_num(Walker* walker, int i, int d) const {
    
    double val = get_val(walker);

    double h = 0.0000001;

    walker->r(i, d) += h;
    walker->make_rel_matrix();
    double val_p = get_val(walker);
    
    walker->r(i, d) -= 2 * h;
    walker->make_rel_matrix();
    double val_m = get_val(walker);

    double deriv = (val_p - val_m)/(2*h*val);
    
    walker->r(i, d) += h;
    walker->make_rel_matrix();
    
    return deriv;
}

double Jastrow::get_laplaciansum_num(Walker* walker) const {

    double val_p, val_m;
    double val = get_val(walker);
    
    double h = 0.00001;
    
    double lapl = 0;

    for (int i = 0; i < n_p; i++){
        for (int j = 0; j < dim; j++){
            
            walker->r(i, j) += h;
            walker->make_rel_matrix();  
            val_p = get_val(walker);
            
            walker->r(i, j) -= 2*h;
            walker->make_rel_matrix();
            val_m = get_val(walker);
            
            walker->r(i, j) += h;
            walker->make_rel_matrix();
            
            lapl += val_p + val_m;
            
        }
    }

    lapl = (lapl - 2*dim*n_p*val)/(h*h*val);
    return lapl;
    
}