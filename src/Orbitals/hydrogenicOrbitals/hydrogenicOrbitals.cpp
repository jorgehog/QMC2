/* 
 * File:   hydrogenicOrbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "hydrogenicOrbitals.h"

#include "../../structs.h"

#include "../../Walker/Walker.h"
#include "../../BasisFunctions/hydrogenic/hydrogenic.h"

#include <iostream>

hydrogenicOrbitals::hydrogenicOrbitals(GeneralParams & gP, VariationalParams & vP)
: Orbitals(gP.n_p, gP.dim) {

    name = "Atoms";

    if (dim != 3) {
        std::cout << "Hydrogenic orbitals are threedimensional creatures." << std::endl;
        exit(1);
    }

    this->alpha = new double();
    this->k = new double();
    this->k2 = new double();

    double* r22d = new double();
    double* r2d = new double();

    this->exp_factor_n1 = new double();
    this->exp_factor_n2 = new double();
    this->exp_factor_n3 = new double();
    this->exp_factor_n4 = new double();

    this->Z = n_p;

    set_parameter(vP.alpha, 0);

    get_qnums();

    basis_functions[0] = new hydrogenic_0(k, k2, r22d, r2d, exp_factor_n1);
    basis_functions[1] = new hydrogenic_1(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[2] = new hydrogenic_2(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[3] = new hydrogenic_3(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[4] = new hydrogenic_4(k, k2, r22d, r2d, exp_factor_n2);
    basis_functions[5] = new hydrogenic_5(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[6] = new hydrogenic_6(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[7] = new hydrogenic_7(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[8] = new hydrogenic_8(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[9] = new hydrogenic_9(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[10] = new hydrogenic_10(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[11] = new hydrogenic_11(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[12] = new hydrogenic_12(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[13] = new hydrogenic_13(k, k2, r22d, r2d, exp_factor_n3);
    basis_functions[14] = new hydrogenic_14(k, k2, r22d, r2d, exp_factor_n4);
    basis_functions[15] = new hydrogenic_15(k, k2, r22d, r2d, exp_factor_n4);
    basis_functions[16] = new hydrogenic_16(k, k2, r22d, r2d, exp_factor_n4);
    basis_functions[17] = new hydrogenic_17(k, k2, r22d, r2d, exp_factor_n4);

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
    dell_basis_functions[0][5] = new dell_hydrogenic_5_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][5] = new dell_hydrogenic_5_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][5] = new dell_hydrogenic_5_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][6] = new dell_hydrogenic_6_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][6] = new dell_hydrogenic_6_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][6] = new dell_hydrogenic_6_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][7] = new dell_hydrogenic_7_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][7] = new dell_hydrogenic_7_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][7] = new dell_hydrogenic_7_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][8] = new dell_hydrogenic_8_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][8] = new dell_hydrogenic_8_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][8] = new dell_hydrogenic_8_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][9] = new dell_hydrogenic_9_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][9] = new dell_hydrogenic_9_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][9] = new dell_hydrogenic_9_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][10] = new dell_hydrogenic_10_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][10] = new dell_hydrogenic_10_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][10] = new dell_hydrogenic_10_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][11] = new dell_hydrogenic_11_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][11] = new dell_hydrogenic_11_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][11] = new dell_hydrogenic_11_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][12] = new dell_hydrogenic_12_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][12] = new dell_hydrogenic_12_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][12] = new dell_hydrogenic_12_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][13] = new dell_hydrogenic_13_x(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[1][13] = new dell_hydrogenic_13_y(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[2][13] = new dell_hydrogenic_13_z(k, k2, r22d, r2d, exp_factor_n3);
    dell_basis_functions[0][14] = new dell_hydrogenic_14_x(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[1][14] = new dell_hydrogenic_14_y(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[2][14] = new dell_hydrogenic_14_z(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[0][15] = new dell_hydrogenic_15_x(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[1][15] = new dell_hydrogenic_15_y(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[2][15] = new dell_hydrogenic_15_z(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[0][16] = new dell_hydrogenic_16_x(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[1][16] = new dell_hydrogenic_16_y(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[2][16] = new dell_hydrogenic_16_z(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[0][17] = new dell_hydrogenic_17_x(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[1][17] = new dell_hydrogenic_17_y(k, k2, r22d, r2d, exp_factor_n4);
    dell_basis_functions[2][17] = new dell_hydrogenic_17_z(k, k2, r22d, r2d, exp_factor_n4);

    lapl_basis_functions[0] = new lapl_hydrogenic_0(k, k2, r22d, r2d, exp_factor_n1);
    lapl_basis_functions[1] = new lapl_hydrogenic_1(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[2] = new lapl_hydrogenic_2(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[3] = new lapl_hydrogenic_3(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[4] = new lapl_hydrogenic_4(k, k2, r22d, r2d, exp_factor_n2);
    lapl_basis_functions[5] = new lapl_hydrogenic_5(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[6] = new lapl_hydrogenic_6(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[7] = new lapl_hydrogenic_7(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[8] = new lapl_hydrogenic_8(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[9] = new lapl_hydrogenic_9(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[10] = new lapl_hydrogenic_10(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[11] = new lapl_hydrogenic_11(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[12] = new lapl_hydrogenic_12(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[13] = new lapl_hydrogenic_13(k, k2, r22d, r2d, exp_factor_n3);
    lapl_basis_functions[14] = new lapl_hydrogenic_14(k, k2, r22d, r2d, exp_factor_n4);
    lapl_basis_functions[15] = new lapl_hydrogenic_15(k, k2, r22d, r2d, exp_factor_n4);
    lapl_basis_functions[16] = new lapl_hydrogenic_16(k, k2, r22d, r2d, exp_factor_n4);
    lapl_basis_functions[17] = new lapl_hydrogenic_17(k, k2, r22d, r2d, exp_factor_n4);

}

void hydrogenicOrbitals::set_qnum_indie_terms(Walker* walker, int i) {

    walker->calc_r_i(i);

    double kr = -(*k) * walker->get_r_i(i);
    *exp_factor_n1 = exp(kr);
    if (nCap > 2) *exp_factor_n2 = exp(kr / 2);
    if (nCap > 10) *exp_factor_n3 = exp(kr / 3);
    if (nCap > 28) *exp_factor_n4 = exp(kr / 4);

}

double hydrogenicOrbitals::get_dell_alpha_phi(Walker* walker, int i, int qnum, int n) {

    (void) n;

    double dphi;

    if (qnum == 0) {

        //-Z*r

        dphi = -Z * walker->get_r_i(i);

    } else if (qnum == 1) {

        //-Z*r*(k*r - 4)/(2*(k*r - 2))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 4) / (2 * ((*k) * walker->get_r_i(i) - 2));

    } else if (qnum == 2) {

        //-Z*r/2

        dphi = -Z * walker->get_r_i(i) / 2;

    } else if (qnum == 3) {

        //-Z*r/2

        dphi = -Z * walker->get_r_i(i) / 2;

    } else if (qnum == 4) {

        //-Z*r/2

        dphi = -Z * walker->get_r_i(i) / 2;

    } else if (qnum == 5) {

        //-Z*r*(2*k^2*r^2 - 30*k*r + 81)/(3*(2*k^2*r^2 - 18*k*r + 27))

        dphi = -Z * walker->get_r_i(i)*(2 * (*k2) * walker->get_r_i2(i) - 30 * (*k) * walker->get_r_i(i) + 81) / (3 * (2 * (*k2) * walker->get_r_i2(i) - 18 * (*k) * walker->get_r_i(i) + 27));

    } else if (qnum == 6) {

        //-Z*r*(k*r - 9)/(3*(k*r - 6))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 9) / (3 * ((*k) * walker->get_r_i(i) - 6));

    } else if (qnum == 7) {

        //-Z*r*(k*r - 9)/(3*(k*r - 6))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 9) / (3 * ((*k) * walker->get_r_i(i) - 6));

    } else if (qnum == 8) {

        //-Z*r*(k*r - 9)/(3*(k*r - 6))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 9) / (3 * ((*k) * walker->get_r_i(i) - 6));

    } else if (qnum == 9) {

        //-Z*r/3

        dphi = -Z * walker->get_r_i(i) / 3;

    } else if (qnum == 10) {

        //-Z*r/3

        dphi = -Z * walker->get_r_i(i) / 3;

    } else if (qnum == 11) {

        //-Z*r/3

        dphi = -Z * walker->get_r_i(i) / 3;

    } else if (qnum == 12) {

        //-Z*r*(x - y)*(x + y)/(3*(x^2 - y^2))

        dphi = -Z * walker->get_r_i(i) / 3;

    } else if (qnum == 13) {

        //-Z*r/3

        dphi = -Z * walker->get_r_i(i) / 3;

    } else if (qnum == 14) {

        //-Z*r*(k^3*r^3 - 36*k^2*r^2 + 336*k*r - 768)/(4*(k^3*r^3 - 24*k^2*r^2 + 144*k*r - 192))

        dphi = -Z * walker->get_r_i(i)*((*k2)*(*k) * walker->get_r_i2(i) * walker->get_r_i(i) - 36 * (*k2) * walker->get_r_i2(i) + 336 * (*k) * walker->get_r_i(i) - 768) / (4 * ((*k2)*(*k) * walker->get_r_i2(i) * walker->get_r_i(i) - 24 * (*k2) * walker->get_r_i2(i) + 144 * (*k) * walker->get_r_i(i) - 192));

    } else if (qnum == 15) {

        //-Z*r*(k*r - 20)*(k*r - 8)/(4*(k^2*r^2 - 20*k*r + 80))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 20)*((*k) * walker->get_r_i(i) - 8) / (4 * ((*k2) * walker->get_r_i2(i) - 20 * (*k) * walker->get_r_i(i) + 80));

    } else if (qnum == 16) {

        //-Z*r*(k*r - 20)*(k*r - 8)/(4*(k^2*r^2 - 20*k*r + 80))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 20)*((*k) * walker->get_r_i(i) - 8) / (4 * ((*k2) * walker->get_r_i2(i) - 20 * (*k) * walker->get_r_i(i) + 80));

    } else if (qnum == 17) {

        //-Z*r*(k*r - 20)*(k*r - 8)/(4*(k^2*r^2 - 20*k*r + 80))

        dphi = -Z * walker->get_r_i(i)*((*k) * walker->get_r_i(i) - 20)*((*k) * walker->get_r_i(i) - 8) / (4 * ((*k2) * walker->get_r_i2(i) - 20 * (*k) * walker->get_r_i(i) + 80));

    } else {
        std::cout << "qnum level " << qnum << " not implemented in dalpha hydro" << std::endl;
        dphi = 0;
    }

    return dphi;

}


double hydrogenicOrbitals::get_sp_energy(int qnum) const {
    int n = qnums(qnum, 0);

    return -Z * Z / (2.0 * n * n);
}

void hydrogenicOrbitals::get_qnums() {
    int n2 = 20;
    qnums.set_size(n2, dim);

    int i = 0;
    int n = 1;
    while (i < n2) {

        for (int l = 0; l < n; l++) {

            if (i == n2) break;

            qnums(i, 0) = n;
            qnums(i, 1) = l;
            qnums(i, 2) = 0;
            i++;

            if (i == n2) break;

            for (int m = 1; m <= l; m++) {
                qnums(i, 0) = n;
                qnums(i, 1) = l;
                qnums(i, 2) = m;
                i++;

                if (i == n2) break;

                qnums(i, 0) = n;
                qnums(i, 1) = l;
                qnums(i, 2) = -m;
                i++;
            }

        }

        n++;
    }

}



