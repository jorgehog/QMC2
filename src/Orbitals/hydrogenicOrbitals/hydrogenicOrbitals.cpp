/* 
 * File:   hydrogenicOrbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "../../QMCheaders.h"

hydrogenicOrbitals::hydrogenicOrbitals(GeneralParams& gp) {
    VariationalParams dummyVar;
    dummyVar.alpha = 1;
    hydrogenicOrbitals(gp, dummyVar);
}

hydrogenicOrbitals::hydrogenicOrbitals(GeneralParams & gP, VariationalParams & vP)
: Orbitals(gP.n_p, gP.dim) {

    this->alpha = new double();
    this->k = new double();
    this->k2 = new double();
    
    this->r22d = new double();
    this->r2d = new double();
    
    this->exp_factor_n1 = new double();
    this->exp_factor_n2 = new double();

    this->Z = (int) gP.systemConstant;
    set_parameter(vP.alpha, 0);
    get_qnums();

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
    *r22d = walker->get_r_i2(i) - walker->r(i, 2) * walker->r(i, 2);
    *r2d = sqrt((*r22d));

}

void hydrogenicOrbitals::get_qnums() {
    //    qnums = arma::zeros<arma::Mat<int> > (n2, dim);
    //    double n_x, n_y;
    //
    //    int n_shells = (int) (0.5 * (sqrt(1 + 4 * n_p) - 1));
    //
    //    int q = 0;
    //    for (int shell = 0; shell < n_shells; shell++) {
    //
    //        n_x = 0;
    //        n_y = shell;
    //
    //        for (int subshell_i = 0; subshell_i <= shell; subshell_i++) {
    //            qnums(q, 0) = n_x;
    //            qnums(q, 1) = n_y;
    //
    //            n_x++;
    //            n_y--;
    //
    //            q++;
    //        }
    //    }
}

double hydrogenicOrbitals::L(int n, int l, double r) const {
    //    if (n < 0) {
    //        return 0;
    //    } else if (n == 0) {
    //        return 1;
    //    } else if (n == 1) {
    //        return 2 * x;
    //    } else if (n == 2) {
    //        return 4 * x * x - 2;
    //    } else if (n == 3) {
    //        return 8 * x * x * x - 12 * x;
    //    } else if (n == 4) {
    //        double x2 = x*x;
    //        return 16 * x2 * x2 - 48 * x2 + 12;
    //    } else {
    //        std::cout << "Unsopported Hermite polynomial level: " << n << std::endl;
    //        exit(1);
    //    }
}

double hydrogenicOrbitals::get_variational_derivative(const Walker* walker, int n) {
    //    double dalpha, sq_w_over_a, exp_fac, H_fac, rij;
    //    int nij;
    //
    //    dalpha = 0;
    //    sq_w_over_a = (*k) / (*alpha);
    //
    //    for (int i = 0; i < n_p; i++) {
    //        exp_fac = -0.5 * w * walker->get_r_i2(i);
    //
    //        for (int qnum = 0; qnum < n2; qnum++) {
    //
    //            H_fac = 0;
    //            for (int j = 0; j < dim; j++) {
    //                rij = walker->r(i, j);
    //                nij = qnums(qnum, j);
    //
    //                H_fac += rij * nij * H(nij - 1, (*k) * rij) / H(nij, (*k) * rij);
    //            }
    //            H_fac *= sq_w_over_a;
    //
    //
    //            dalpha += walker->inv(qnum, i)*(H_fac + exp_fac) * walker->phi(i, qnum);
    //        }
    //    }
    //
    //    return dalpha;
    return 0;
}