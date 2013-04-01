/* 
 * File:   hydrogenicOrbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "../../QMCheaders.h"

hydrogenicOrbitals::hydrogenicOrbitals(GeneralParams & gP, VariationalParams & vP)
: Orbitals(gP.n_p, gP.dim) {

    this->alpha = new double();
    this->k = new double();
    this->k2 = new double();

    double* r22d = new double();
    double* r2d = new double();

    this->exp_factor_n1 = new double();
    this->exp_factor_n2 = new double();

    this->Z = (int) gP.systemConstant;
    set_parameter(vP.alpha, 0);

    basis_functions[0] = new hydrogenic_0(k, k2, r22d, r2d, exp_factor_n1);
    basis_functions[1] = new hydrogenic_1(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[2] = new hydrogenic_2(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[3] = new hydrogenic_3(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[4] = new hydrogenic_4(k, k2, r22d, r2d, exp_factor_n2);

    dell_basis_functions[0][0] = new dell_hydrogenic_0_x(k, k2, r22d, r2d, exp_factor_n1);
    dell_basis_functions[1][0] = new dell_hydrogenic_0_y(k, k2, r22d, r2d, exp_factor_n1);
    dell_basis_functions[2][0] = new dell_hydrogenic_0_z(k, k2, r22d, r2d, exp_factor_n1);
    dell_basis_functions[0][1] = new dell_hydrogenic_1_x(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[1][1] = new dell_hydrogenic_1_y(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[2][1] = new dell_hydrogenic_1_z(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[0][2] = new dell_hydrogenic_2_x(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[1][2] = new dell_hydrogenic_2_y(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[2][2] = new dell_hydrogenic_2_z(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[0][3] = new dell_hydrogenic_3_x(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[1][3] = new dell_hydrogenic_3_y(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[2][3] = new dell_hydrogenic_3_z(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[0][4] = new dell_hydrogenic_4_x(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[1][4] = new dell_hydrogenic_4_y(k, k2, r22d, r2d, exp_factor_n2);
    dell_basis_functions[2][4] = new dell_hydrogenic_4_z(k, k2, r22d, r2d, exp_factor_n2);

    lapl_basis_functions[0] = new lapl_hydrogenic_0(k, k2, r22d, r2d, exp_factor_n1);
    lapl_basis_functions[1] = new lapl_hydrogenic_1(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[2] = new lapl_hydrogenic_2(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[3] = new lapl_hydrogenic_3(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[4] = new lapl_hydrogenic_4(k, k2, r22d, r2d, exp_factor_n2);
}

void hydrogenicOrbitals::set_qnum_indie_terms(const Walker* walker, int i) {

    *exp_factor_n1 = exp(-(*k) * walker->get_r_i(i));
    *exp_factor_n2 = exp(-(*k) * walker->get_r_i(i) / 2);

}
double hydrogenicOrbitals::get_dell_alpha_phi(const Walker* walker, int qnum, int i){
    
    double dphi;
    
    if (qnum == 0) {    
    
        //-Z*r
        
        dphi = -Z*walker->get_r_i(i);
        
    } else if (qnum == 1) {    
    
        //-Z*r*(k*r - 4)/(2*(k*r - 2))
        
        dphi = -Z*walker->get_r_i(i)*((*k)*walker->get_r_i(i) - 4)/(2*((*k)*walker->get_r_i(i) - 2));
        
    } else if (qnum == 2) {    
    
        //-Z*r/2
        
        dphi = -Z*walker->get_r_i(i)/2;
        
    } else if (qnum == 3) {    
    
        //-Z*r/2
        
        dphi = -Z*walker->get_r_i(i)/2;
        
    } else if (qnum == 4) {    
    
        //-Z*r/2
        
        dphi = -Z*walker->get_r_i(i)/2;
        
    }
    
    return dphi;
    
}


double hydrogenicOrbitals::get_variational_derivative(const Walker* walker, int n) {
    double dalpha, dell_alpha_phi;

    dalpha = 0;

    for (int i = 0; i < n_p; i++) {
        for (int qnum = 0; qnum < n2; qnum++) {

            dell_alpha_phi = get_dell_alpha_phi(walker, qnum, i);
            dalpha += walker->inv(qnum, i) * dell_alpha_phi * walker->phi(i, qnum);

        }
    }

    return dalpha;

}
