/* 
 * File:   AlphaHarmonicOscillator.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "../../QMCheaders.h"

AlphaHarmonicOscillator::AlphaHarmonicOscillator(GeneralParams& gp) {
    VariationalParams dummyVar;
    dummyVar.alpha = 1;
    AlphaHarmonicOscillator(gp, dummyVar);
}

AlphaHarmonicOscillator::AlphaHarmonicOscillator(GeneralParams & gP, VariationalParams & vP)
: Orbitals(gP.n_p, gP.dim) {

    this->alpha = new double();
    this->k = new double();
    this->k2 = new double();

    this->w = gP.w;
    set_parameter(vP.alpha, 0);
    get_qnums();

    basis_functions[0] = new alphaHO_0(k, k2);
    basis_functions[1] = new alphaHO_1(k, k2);
    basis_functions[2] = new alphaHO_2(k, k2);
    basis_functions[3] = new alphaHO_3(k, k2);
    basis_functions[4] = new alphaHO_4(k, k2);
    basis_functions[5] = new alphaHO_5(k, k2);
    basis_functions[6] = new alphaHO_6(k, k2);
    basis_functions[7] = new alphaHO_7(k, k2);
    basis_functions[8] = new alphaHO_8(k, k2);
    basis_functions[9] = new alphaHO_9(k, k2);
    basis_functions[10] = new alphaHO_10(k, k2);
    basis_functions[11] = new alphaHO_11(k, k2);
    basis_functions[12] = new alphaHO_12(k, k2);
    basis_functions[13] = new alphaHO_13(k, k2);
    basis_functions[14] = new alphaHO_14(k, k2);

    dell_basis_functions[0][0] = new dell_alphaHO_0_x(k, k2);
    dell_basis_functions[1][0] = new dell_alphaHO_0_y(k, k2);
    dell_basis_functions[0][1] = new dell_alphaHO_1_x(k, k2);
    dell_basis_functions[1][1] = new dell_alphaHO_1_y(k, k2);
    dell_basis_functions[0][2] = new dell_alphaHO_2_x(k, k2);
    dell_basis_functions[1][2] = new dell_alphaHO_2_y(k, k2);
    dell_basis_functions[0][3] = new dell_alphaHO_3_x(k, k2);
    dell_basis_functions[1][3] = new dell_alphaHO_3_y(k, k2);
    dell_basis_functions[0][4] = new dell_alphaHO_4_x(k, k2);
    dell_basis_functions[1][4] = new dell_alphaHO_4_y(k, k2);
    dell_basis_functions[0][5] = new dell_alphaHO_5_x(k, k2);
    dell_basis_functions[1][5] = new dell_alphaHO_5_y(k, k2);
    dell_basis_functions[0][6] = new dell_alphaHO_6_x(k, k2);
    dell_basis_functions[1][6] = new dell_alphaHO_6_y(k, k2);
    dell_basis_functions[0][7] = new dell_alphaHO_7_x(k, k2);
    dell_basis_functions[1][7] = new dell_alphaHO_7_y(k, k2);
    dell_basis_functions[0][8] = new dell_alphaHO_8_x(k, k2);
    dell_basis_functions[1][8] = new dell_alphaHO_8_y(k, k2);
    dell_basis_functions[0][9] = new dell_alphaHO_9_x(k, k2);
    dell_basis_functions[1][9] = new dell_alphaHO_9_y(k, k2);
    dell_basis_functions[0][10] = new dell_alphaHO_10_x(k, k2);
    dell_basis_functions[1][10] = new dell_alphaHO_10_y(k, k2);
    dell_basis_functions[0][11] = new dell_alphaHO_11_x(k, k2);
    dell_basis_functions[1][11] = new dell_alphaHO_11_y(k, k2);
    dell_basis_functions[0][12] = new dell_alphaHO_12_x(k, k2);
    dell_basis_functions[1][12] = new dell_alphaHO_12_y(k, k2);
    dell_basis_functions[0][13] = new dell_alphaHO_13_x(k, k2);
    dell_basis_functions[1][13] = new dell_alphaHO_13_y(k, k2);
    dell_basis_functions[0][14] = new dell_alphaHO_14_x(k, k2);
    dell_basis_functions[1][14] = new dell_alphaHO_14_y(k, k2);


    lapl_basis_functions[0] = new lapl_alphaHO_0(k, k2);
    lapl_basis_functions[1] = new lapl_alphaHO_1(k, k2);
    lapl_basis_functions[2] = new lapl_alphaHO_2(k, k2);
    lapl_basis_functions[3] = new lapl_alphaHO_3(k, k2);
    lapl_basis_functions[4] = new lapl_alphaHO_4(k, k2);
    lapl_basis_functions[5] = new lapl_alphaHO_5(k, k2);
    lapl_basis_functions[6] = new lapl_alphaHO_6(k, k2);
    lapl_basis_functions[7] = new lapl_alphaHO_7(k, k2);
    lapl_basis_functions[8] = new lapl_alphaHO_8(k, k2);
    lapl_basis_functions[9] = new lapl_alphaHO_9(k, k2);
    lapl_basis_functions[10] = new lapl_alphaHO_10(k, k2);
    lapl_basis_functions[11] = new lapl_alphaHO_11(k, k2);
    lapl_basis_functions[12] = new lapl_alphaHO_12(k, k2);
    lapl_basis_functions[13] = new lapl_alphaHO_13(k, k2);
    lapl_basis_functions[14] = new lapl_alphaHO_14(k, k2);
}

double AlphaHarmonicOscillator::get_parameter(int n) {
    return *alpha;
}

void AlphaHarmonicOscillator::set_parameter(double parameter, int n) {
    *alpha = parameter;
    *k2 = parameter*w;
    *k = sqrt(*k2);
}

void AlphaHarmonicOscillator::get_qnums() {
    qnums = arma::zeros<arma::Mat<int> > (n2, dim);
    double n_x, n_y;

    int n_shells = (int) (0.5 * (sqrt(1 + 4 * n_p) - 1));

    int q = 0;
    for (int shell = 0; shell < n_shells; shell++) {

        n_x = 0;
        n_y = shell;

        for (int subshell_i = 0; subshell_i <= shell; subshell_i++) {
            qnums(q, 0) = n_x;
            qnums(q, 1) = n_y;

            n_x++;
            n_y--;

            q++;
        }
    }
}

double AlphaHarmonicOscillator::H(int n, double x) const {
    if (n < 0) {
        return 0;
    } else if (n == 0) {
        return 1;
    } else if (n == 1) {
        return 2 * x;
    } else if (n == 2) {
        return 4 * x * x - 2;
    } else if (n == 3) {
        return 8 * x * x * x - 12 * x;
    } else if (n == 4) {
        double x2 = x*x;
        return 16 * x2 * x2 - 48 * x2 + 12;
    } else {
        std::cout << "Unsopported Hermite polynomial level: " << n << std::endl;
        exit(1);
    }
}

double AlphaHarmonicOscillator::get_variational_derivative(const Walker* walker, int n) const {
    double dalpha, sq_w_over_a, exp_fac, H_fac, rij;
    int nij;

    dalpha = 0;
    sq_w_over_a = (*k) / (*alpha);

    for (int i = 0; i < n_p; i++) {
        exp_fac = -0.5 * w * walker->get_r_i2(i);

        for (int qnum = 0; qnum < n2; qnum++) {

            H_fac = 0;
            for (int j = 0; j < dim; j++) {
                rij = walker->r(i, j);
                nij = qnums(qnum, j);

                H_fac += rij * nij * H(nij - 1, (*k) * rij) / H(nij, (*k) * rij);
            }
            H_fac *= sq_w_over_a;


            dalpha += walker->inv(i, qnum)*(H_fac + exp_fac) * basis_functions[qnum]->eval(walker, i);
        }
    }


    return dalpha;
}