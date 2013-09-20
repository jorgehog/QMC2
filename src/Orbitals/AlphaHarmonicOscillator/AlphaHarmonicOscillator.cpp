/* 
 * File:   AlphaHarmonicOscillator.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "AlphaHarmonicOscillator.h"

#include "../../structs.h"

#include "../../Walker/Walker.h"

#include "../../BasisFunctions/HarmonicOscillator/HarmonicOscillator.h"
#include "../../BasisFunctions/HarmonicOscillator3D/HarmonicOscillator3D.h"


AlphaHarmonicOscillator::AlphaHarmonicOscillator(GeneralParams & gP, VariationalParams & vP)
: Orbitals(gP.n_p, gP.dim) {

    this->alpha = new double();
    this->k = new double();
    this->k2 = new double();
    this->exp_factor = new double();

    this->w = gP.systemConstant;
    set_parameter(vP.alpha, 0);

    qnums = arma::zeros<arma::imat > (n2, dim);

    if (dim == 2) {
        name = "QDots";
        get_qnums();
        setup_basis();
    } else if (dim == 3) {
        name = "QDots3D";
        get_qnums3D();
        setup_basis3D();
    }
}

void AlphaHarmonicOscillator::setup_basis3D() {

    basis_functions[0] = new HarmonicOscillator3D_0(k, k2, exp_factor);
    basis_functions[1] = new HarmonicOscillator3D_1(k, k2, exp_factor);
    basis_functions[2] = new HarmonicOscillator3D_2(k, k2, exp_factor);
    basis_functions[3] = new HarmonicOscillator3D_3(k, k2, exp_factor);
    basis_functions[4] = new HarmonicOscillator3D_4(k, k2, exp_factor);
    basis_functions[5] = new HarmonicOscillator3D_5(k, k2, exp_factor);
    basis_functions[6] = new HarmonicOscillator3D_6(k, k2, exp_factor);
    basis_functions[7] = new HarmonicOscillator3D_7(k, k2, exp_factor);
    basis_functions[8] = new HarmonicOscillator3D_8(k, k2, exp_factor);
    basis_functions[9] = new HarmonicOscillator3D_9(k, k2, exp_factor);
    basis_functions[10] = new HarmonicOscillator3D_10(k, k2, exp_factor);
    basis_functions[11] = new HarmonicOscillator3D_11(k, k2, exp_factor);
    basis_functions[12] = new HarmonicOscillator3D_12(k, k2, exp_factor);
    basis_functions[13] = new HarmonicOscillator3D_13(k, k2, exp_factor);
    basis_functions[14] = new HarmonicOscillator3D_14(k, k2, exp_factor);
    basis_functions[15] = new HarmonicOscillator3D_15(k, k2, exp_factor);
    basis_functions[16] = new HarmonicOscillator3D_16(k, k2, exp_factor);
    basis_functions[17] = new HarmonicOscillator3D_17(k, k2, exp_factor);
    basis_functions[18] = new HarmonicOscillator3D_18(k, k2, exp_factor);
    basis_functions[19] = new HarmonicOscillator3D_19(k, k2, exp_factor);

    dell_basis_functions[0][0] = new dell_HarmonicOscillator3D_0_x(k, k2, exp_factor);
    dell_basis_functions[1][0] = new dell_HarmonicOscillator3D_0_y(k, k2, exp_factor);
    dell_basis_functions[2][0] = new dell_HarmonicOscillator3D_0_z(k, k2, exp_factor);
    dell_basis_functions[0][1] = new dell_HarmonicOscillator3D_1_x(k, k2, exp_factor);
    dell_basis_functions[1][1] = new dell_HarmonicOscillator3D_1_y(k, k2, exp_factor);
    dell_basis_functions[2][1] = new dell_HarmonicOscillator3D_1_z(k, k2, exp_factor);
    dell_basis_functions[0][2] = new dell_HarmonicOscillator3D_2_x(k, k2, exp_factor);
    dell_basis_functions[1][2] = new dell_HarmonicOscillator3D_2_y(k, k2, exp_factor);
    dell_basis_functions[2][2] = new dell_HarmonicOscillator3D_2_z(k, k2, exp_factor);
    dell_basis_functions[0][3] = new dell_HarmonicOscillator3D_3_x(k, k2, exp_factor);
    dell_basis_functions[1][3] = new dell_HarmonicOscillator3D_3_y(k, k2, exp_factor);
    dell_basis_functions[2][3] = new dell_HarmonicOscillator3D_3_z(k, k2, exp_factor);
    dell_basis_functions[0][4] = new dell_HarmonicOscillator3D_4_x(k, k2, exp_factor);
    dell_basis_functions[1][4] = new dell_HarmonicOscillator3D_4_y(k, k2, exp_factor);
    dell_basis_functions[2][4] = new dell_HarmonicOscillator3D_4_z(k, k2, exp_factor);
    dell_basis_functions[0][5] = new dell_HarmonicOscillator3D_5_x(k, k2, exp_factor);
    dell_basis_functions[1][5] = new dell_HarmonicOscillator3D_5_y(k, k2, exp_factor);
    dell_basis_functions[2][5] = new dell_HarmonicOscillator3D_5_z(k, k2, exp_factor);
    dell_basis_functions[0][6] = new dell_HarmonicOscillator3D_6_x(k, k2, exp_factor);
    dell_basis_functions[1][6] = new dell_HarmonicOscillator3D_6_y(k, k2, exp_factor);
    dell_basis_functions[2][6] = new dell_HarmonicOscillator3D_6_z(k, k2, exp_factor);
    dell_basis_functions[0][7] = new dell_HarmonicOscillator3D_7_x(k, k2, exp_factor);
    dell_basis_functions[1][7] = new dell_HarmonicOscillator3D_7_y(k, k2, exp_factor);
    dell_basis_functions[2][7] = new dell_HarmonicOscillator3D_7_z(k, k2, exp_factor);
    dell_basis_functions[0][8] = new dell_HarmonicOscillator3D_8_x(k, k2, exp_factor);
    dell_basis_functions[1][8] = new dell_HarmonicOscillator3D_8_y(k, k2, exp_factor);
    dell_basis_functions[2][8] = new dell_HarmonicOscillator3D_8_z(k, k2, exp_factor);
    dell_basis_functions[0][9] = new dell_HarmonicOscillator3D_9_x(k, k2, exp_factor);
    dell_basis_functions[1][9] = new dell_HarmonicOscillator3D_9_y(k, k2, exp_factor);
    dell_basis_functions[2][9] = new dell_HarmonicOscillator3D_9_z(k, k2, exp_factor);
    dell_basis_functions[0][10] = new dell_HarmonicOscillator3D_10_x(k, k2, exp_factor);
    dell_basis_functions[1][10] = new dell_HarmonicOscillator3D_10_y(k, k2, exp_factor);
    dell_basis_functions[2][10] = new dell_HarmonicOscillator3D_10_z(k, k2, exp_factor);
    dell_basis_functions[0][11] = new dell_HarmonicOscillator3D_11_x(k, k2, exp_factor);
    dell_basis_functions[1][11] = new dell_HarmonicOscillator3D_11_y(k, k2, exp_factor);
    dell_basis_functions[2][11] = new dell_HarmonicOscillator3D_11_z(k, k2, exp_factor);
    dell_basis_functions[0][12] = new dell_HarmonicOscillator3D_12_x(k, k2, exp_factor);
    dell_basis_functions[1][12] = new dell_HarmonicOscillator3D_12_y(k, k2, exp_factor);
    dell_basis_functions[2][12] = new dell_HarmonicOscillator3D_12_z(k, k2, exp_factor);
    dell_basis_functions[0][13] = new dell_HarmonicOscillator3D_13_x(k, k2, exp_factor);
    dell_basis_functions[1][13] = new dell_HarmonicOscillator3D_13_y(k, k2, exp_factor);
    dell_basis_functions[2][13] = new dell_HarmonicOscillator3D_13_z(k, k2, exp_factor);
    dell_basis_functions[0][14] = new dell_HarmonicOscillator3D_14_x(k, k2, exp_factor);
    dell_basis_functions[1][14] = new dell_HarmonicOscillator3D_14_y(k, k2, exp_factor);
    dell_basis_functions[2][14] = new dell_HarmonicOscillator3D_14_z(k, k2, exp_factor);
    dell_basis_functions[0][15] = new dell_HarmonicOscillator3D_15_x(k, k2, exp_factor);
    dell_basis_functions[1][15] = new dell_HarmonicOscillator3D_15_y(k, k2, exp_factor);
    dell_basis_functions[2][15] = new dell_HarmonicOscillator3D_15_z(k, k2, exp_factor);
    dell_basis_functions[0][16] = new dell_HarmonicOscillator3D_16_x(k, k2, exp_factor);
    dell_basis_functions[1][16] = new dell_HarmonicOscillator3D_16_y(k, k2, exp_factor);
    dell_basis_functions[2][16] = new dell_HarmonicOscillator3D_16_z(k, k2, exp_factor);
    dell_basis_functions[0][17] = new dell_HarmonicOscillator3D_17_x(k, k2, exp_factor);
    dell_basis_functions[1][17] = new dell_HarmonicOscillator3D_17_y(k, k2, exp_factor);
    dell_basis_functions[2][17] = new dell_HarmonicOscillator3D_17_z(k, k2, exp_factor);
    dell_basis_functions[0][18] = new dell_HarmonicOscillator3D_18_x(k, k2, exp_factor);
    dell_basis_functions[1][18] = new dell_HarmonicOscillator3D_18_y(k, k2, exp_factor);
    dell_basis_functions[2][18] = new dell_HarmonicOscillator3D_18_z(k, k2, exp_factor);
    dell_basis_functions[0][19] = new dell_HarmonicOscillator3D_19_x(k, k2, exp_factor);
    dell_basis_functions[1][19] = new dell_HarmonicOscillator3D_19_y(k, k2, exp_factor);
    dell_basis_functions[2][19] = new dell_HarmonicOscillator3D_19_z(k, k2, exp_factor);

    lapl_basis_functions[0] = new lapl_HarmonicOscillator3D_0(k, k2, exp_factor);
    lapl_basis_functions[1] = new lapl_HarmonicOscillator3D_1(k, k2, exp_factor);
    lapl_basis_functions[2] = new lapl_HarmonicOscillator3D_2(k, k2, exp_factor);
    lapl_basis_functions[3] = new lapl_HarmonicOscillator3D_3(k, k2, exp_factor);
    lapl_basis_functions[4] = new lapl_HarmonicOscillator3D_4(k, k2, exp_factor);
    lapl_basis_functions[5] = new lapl_HarmonicOscillator3D_5(k, k2, exp_factor);
    lapl_basis_functions[6] = new lapl_HarmonicOscillator3D_6(k, k2, exp_factor);
    lapl_basis_functions[7] = new lapl_HarmonicOscillator3D_7(k, k2, exp_factor);
    lapl_basis_functions[8] = new lapl_HarmonicOscillator3D_8(k, k2, exp_factor);
    lapl_basis_functions[9] = new lapl_HarmonicOscillator3D_9(k, k2, exp_factor);
    lapl_basis_functions[10] = new lapl_HarmonicOscillator3D_10(k, k2, exp_factor);
    lapl_basis_functions[11] = new lapl_HarmonicOscillator3D_11(k, k2, exp_factor);
    lapl_basis_functions[12] = new lapl_HarmonicOscillator3D_12(k, k2, exp_factor);
    lapl_basis_functions[13] = new lapl_HarmonicOscillator3D_13(k, k2, exp_factor);
    lapl_basis_functions[14] = new lapl_HarmonicOscillator3D_14(k, k2, exp_factor);
    lapl_basis_functions[15] = new lapl_HarmonicOscillator3D_15(k, k2, exp_factor);
    lapl_basis_functions[16] = new lapl_HarmonicOscillator3D_16(k, k2, exp_factor);
    lapl_basis_functions[17] = new lapl_HarmonicOscillator3D_17(k, k2, exp_factor);
    lapl_basis_functions[18] = new lapl_HarmonicOscillator3D_18(k, k2, exp_factor);
    lapl_basis_functions[19] = new lapl_HarmonicOscillator3D_19(k, k2, exp_factor);

}

void AlphaHarmonicOscillator::get_qnums3D() {

    qnums(0, 0) = 0;
    qnums(0, 1) = 0;
    qnums(0, 2) = 0;

    if (n_p == 2) return;

    qnums(1, 0) = 0;
    qnums(1, 1) = 0;
    qnums(1, 2) = 1;
    qnums(2, 0) = 0;
    qnums(2, 1) = 1;
    qnums(2, 2) = 0;
    qnums(3, 0) = 1;
    qnums(3, 1) = 0;
    qnums(3, 2) = 0;

    if (n_p == 8) return;

    qnums(4, 0) = 0;
    qnums(4, 1) = 0;
    qnums(4, 2) = 2;
    qnums(5, 0) = 0;
    qnums(5, 1) = 1;
    qnums(5, 2) = 1;
    qnums(6, 0) = 0;
    qnums(6, 1) = 2;
    qnums(6, 2) = 0;
    qnums(7, 0) = 1;
    qnums(7, 1) = 0;
    qnums(7, 2) = 1;
    qnums(8, 0) = 1;
    qnums(8, 1) = 1;
    qnums(8, 2) = 0;
    qnums(9, 0) = 2;
    qnums(9, 1) = 0;
    qnums(9, 2) = 0;

    if (n_p == 20) return;

    qnums(10, 0) = 0;
    qnums(10, 1) = 0;
    qnums(10, 2) = 3;
    qnums(11, 0) = 0;
    qnums(11, 1) = 1;
    qnums(11, 2) = 2;
    qnums(12, 0) = 0;
    qnums(12, 1) = 2;
    qnums(12, 2) = 1;
    qnums(13, 0) = 0;
    qnums(13, 1) = 3;
    qnums(13, 2) = 0;
    qnums(14, 0) = 1;
    qnums(14, 1) = 0;
    qnums(14, 2) = 2;
    qnums(15, 0) = 1;
    qnums(15, 1) = 1;
    qnums(15, 2) = 1;
    qnums(16, 0) = 1;
    qnums(16, 1) = 2;
    qnums(16, 2) = 0;
    qnums(17, 0) = 2;
    qnums(17, 1) = 0;
    qnums(17, 2) = 1;
    qnums(18, 0) = 2;
    qnums(18, 1) = 1;
    qnums(18, 2) = 0;
    qnums(19, 0) = 3;
    qnums(19, 1) = 0;
    qnums(19, 2) = 0;

    if (n_p != 40) {
        std::cout << "Quantum dots open shell symmetry broken or N > 40" << std::endl;
        exit(1);
    }

}

void AlphaHarmonicOscillator::setup_basis() {

    basis_functions[0] = new HarmonicOscillator_0(k, k2, exp_factor);
    basis_functions[1] = new HarmonicOscillator_1(k, k2, exp_factor);
    basis_functions[2] = new HarmonicOscillator_2(k, k2, exp_factor);
    basis_functions[3] = new HarmonicOscillator_3(k, k2, exp_factor);
    basis_functions[4] = new HarmonicOscillator_4(k, k2, exp_factor);
    basis_functions[5] = new HarmonicOscillator_5(k, k2, exp_factor);
    basis_functions[6] = new HarmonicOscillator_6(k, k2, exp_factor);
    basis_functions[7] = new HarmonicOscillator_7(k, k2, exp_factor);
    basis_functions[8] = new HarmonicOscillator_8(k, k2, exp_factor);
    basis_functions[9] = new HarmonicOscillator_9(k, k2, exp_factor);
    basis_functions[10] = new HarmonicOscillator_10(k, k2, exp_factor);
    basis_functions[11] = new HarmonicOscillator_11(k, k2, exp_factor);
    basis_functions[12] = new HarmonicOscillator_12(k, k2, exp_factor);
    basis_functions[13] = new HarmonicOscillator_13(k, k2, exp_factor);
    basis_functions[14] = new HarmonicOscillator_14(k, k2, exp_factor);
    basis_functions[15] = new HarmonicOscillator_15(k, k2, exp_factor);
    basis_functions[16] = new HarmonicOscillator_16(k, k2, exp_factor);
    basis_functions[17] = new HarmonicOscillator_17(k, k2, exp_factor);
    basis_functions[18] = new HarmonicOscillator_18(k, k2, exp_factor);
    basis_functions[19] = new HarmonicOscillator_19(k, k2, exp_factor);
    basis_functions[20] = new HarmonicOscillator_20(k, k2, exp_factor);
    basis_functions[21] = new HarmonicOscillator_21(k, k2, exp_factor);
    basis_functions[22] = new HarmonicOscillator_22(k, k2, exp_factor);
    basis_functions[23] = new HarmonicOscillator_23(k, k2, exp_factor);
    basis_functions[24] = new HarmonicOscillator_24(k, k2, exp_factor);
    basis_functions[25] = new HarmonicOscillator_25(k, k2, exp_factor);
    basis_functions[26] = new HarmonicOscillator_26(k, k2, exp_factor);
    basis_functions[27] = new HarmonicOscillator_27(k, k2, exp_factor);

    dell_basis_functions[0][0] = new dell_HarmonicOscillator_0_x(k, k2, exp_factor);
    dell_basis_functions[1][0] = new dell_HarmonicOscillator_0_y(k, k2, exp_factor);
    dell_basis_functions[0][1] = new dell_HarmonicOscillator_1_x(k, k2, exp_factor);
    dell_basis_functions[1][1] = new dell_HarmonicOscillator_1_y(k, k2, exp_factor);
    dell_basis_functions[0][2] = new dell_HarmonicOscillator_2_x(k, k2, exp_factor);
    dell_basis_functions[1][2] = new dell_HarmonicOscillator_2_y(k, k2, exp_factor);
    dell_basis_functions[0][3] = new dell_HarmonicOscillator_3_x(k, k2, exp_factor);
    dell_basis_functions[1][3] = new dell_HarmonicOscillator_3_y(k, k2, exp_factor);
    dell_basis_functions[0][4] = new dell_HarmonicOscillator_4_x(k, k2, exp_factor);
    dell_basis_functions[1][4] = new dell_HarmonicOscillator_4_y(k, k2, exp_factor);
    dell_basis_functions[0][5] = new dell_HarmonicOscillator_5_x(k, k2, exp_factor);
    dell_basis_functions[1][5] = new dell_HarmonicOscillator_5_y(k, k2, exp_factor);
    dell_basis_functions[0][6] = new dell_HarmonicOscillator_6_x(k, k2, exp_factor);
    dell_basis_functions[1][6] = new dell_HarmonicOscillator_6_y(k, k2, exp_factor);
    dell_basis_functions[0][7] = new dell_HarmonicOscillator_7_x(k, k2, exp_factor);
    dell_basis_functions[1][7] = new dell_HarmonicOscillator_7_y(k, k2, exp_factor);
    dell_basis_functions[0][8] = new dell_HarmonicOscillator_8_x(k, k2, exp_factor);
    dell_basis_functions[1][8] = new dell_HarmonicOscillator_8_y(k, k2, exp_factor);
    dell_basis_functions[0][9] = new dell_HarmonicOscillator_9_x(k, k2, exp_factor);
    dell_basis_functions[1][9] = new dell_HarmonicOscillator_9_y(k, k2, exp_factor);
    dell_basis_functions[0][10] = new dell_HarmonicOscillator_10_x(k, k2, exp_factor);
    dell_basis_functions[1][10] = new dell_HarmonicOscillator_10_y(k, k2, exp_factor);
    dell_basis_functions[0][11] = new dell_HarmonicOscillator_11_x(k, k2, exp_factor);
    dell_basis_functions[1][11] = new dell_HarmonicOscillator_11_y(k, k2, exp_factor);
    dell_basis_functions[0][12] = new dell_HarmonicOscillator_12_x(k, k2, exp_factor);
    dell_basis_functions[1][12] = new dell_HarmonicOscillator_12_y(k, k2, exp_factor);
    dell_basis_functions[0][13] = new dell_HarmonicOscillator_13_x(k, k2, exp_factor);
    dell_basis_functions[1][13] = new dell_HarmonicOscillator_13_y(k, k2, exp_factor);
    dell_basis_functions[0][14] = new dell_HarmonicOscillator_14_x(k, k2, exp_factor);
    dell_basis_functions[1][14] = new dell_HarmonicOscillator_14_y(k, k2, exp_factor);
    dell_basis_functions[0][15] = new dell_HarmonicOscillator_15_x(k, k2, exp_factor);
    dell_basis_functions[1][15] = new dell_HarmonicOscillator_15_y(k, k2, exp_factor);
    dell_basis_functions[0][16] = new dell_HarmonicOscillator_16_x(k, k2, exp_factor);
    dell_basis_functions[1][16] = new dell_HarmonicOscillator_16_y(k, k2, exp_factor);
    dell_basis_functions[0][17] = new dell_HarmonicOscillator_17_x(k, k2, exp_factor);
    dell_basis_functions[1][17] = new dell_HarmonicOscillator_17_y(k, k2, exp_factor);
    dell_basis_functions[0][18] = new dell_HarmonicOscillator_18_x(k, k2, exp_factor);
    dell_basis_functions[1][18] = new dell_HarmonicOscillator_18_y(k, k2, exp_factor);
    dell_basis_functions[0][19] = new dell_HarmonicOscillator_19_x(k, k2, exp_factor);
    dell_basis_functions[1][19] = new dell_HarmonicOscillator_19_y(k, k2, exp_factor);
    dell_basis_functions[0][20] = new dell_HarmonicOscillator_20_x(k, k2, exp_factor);
    dell_basis_functions[1][20] = new dell_HarmonicOscillator_20_y(k, k2, exp_factor);
    dell_basis_functions[0][21] = new dell_HarmonicOscillator_21_x(k, k2, exp_factor);
    dell_basis_functions[1][21] = new dell_HarmonicOscillator_21_y(k, k2, exp_factor);
    dell_basis_functions[0][22] = new dell_HarmonicOscillator_22_x(k, k2, exp_factor);
    dell_basis_functions[1][22] = new dell_HarmonicOscillator_22_y(k, k2, exp_factor);
    dell_basis_functions[0][23] = new dell_HarmonicOscillator_23_x(k, k2, exp_factor);
    dell_basis_functions[1][23] = new dell_HarmonicOscillator_23_y(k, k2, exp_factor);
    dell_basis_functions[0][24] = new dell_HarmonicOscillator_24_x(k, k2, exp_factor);
    dell_basis_functions[1][24] = new dell_HarmonicOscillator_24_y(k, k2, exp_factor);
    dell_basis_functions[0][25] = new dell_HarmonicOscillator_25_x(k, k2, exp_factor);
    dell_basis_functions[1][25] = new dell_HarmonicOscillator_25_y(k, k2, exp_factor);
    dell_basis_functions[0][26] = new dell_HarmonicOscillator_26_x(k, k2, exp_factor);
    dell_basis_functions[1][26] = new dell_HarmonicOscillator_26_y(k, k2, exp_factor);
    dell_basis_functions[0][27] = new dell_HarmonicOscillator_27_x(k, k2, exp_factor);
    dell_basis_functions[1][27] = new dell_HarmonicOscillator_27_y(k, k2, exp_factor);


    lapl_basis_functions[0] = new lapl_HarmonicOscillator_0(k, k2, exp_factor);
    lapl_basis_functions[1] = new lapl_HarmonicOscillator_1(k, k2, exp_factor);
    lapl_basis_functions[2] = new lapl_HarmonicOscillator_2(k, k2, exp_factor);
    lapl_basis_functions[3] = new lapl_HarmonicOscillator_3(k, k2, exp_factor);
    lapl_basis_functions[4] = new lapl_HarmonicOscillator_4(k, k2, exp_factor);
    lapl_basis_functions[5] = new lapl_HarmonicOscillator_5(k, k2, exp_factor);
    lapl_basis_functions[6] = new lapl_HarmonicOscillator_6(k, k2, exp_factor);
    lapl_basis_functions[7] = new lapl_HarmonicOscillator_7(k, k2, exp_factor);
    lapl_basis_functions[8] = new lapl_HarmonicOscillator_8(k, k2, exp_factor);
    lapl_basis_functions[9] = new lapl_HarmonicOscillator_9(k, k2, exp_factor);
    lapl_basis_functions[10] = new lapl_HarmonicOscillator_10(k, k2, exp_factor);
    lapl_basis_functions[11] = new lapl_HarmonicOscillator_11(k, k2, exp_factor);
    lapl_basis_functions[12] = new lapl_HarmonicOscillator_12(k, k2, exp_factor);
    lapl_basis_functions[13] = new lapl_HarmonicOscillator_13(k, k2, exp_factor);
    lapl_basis_functions[14] = new lapl_HarmonicOscillator_14(k, k2, exp_factor);
    lapl_basis_functions[15] = new lapl_HarmonicOscillator_15(k, k2, exp_factor);
    lapl_basis_functions[16] = new lapl_HarmonicOscillator_16(k, k2, exp_factor);
    lapl_basis_functions[17] = new lapl_HarmonicOscillator_17(k, k2, exp_factor);
    lapl_basis_functions[18] = new lapl_HarmonicOscillator_18(k, k2, exp_factor);
    lapl_basis_functions[19] = new lapl_HarmonicOscillator_19(k, k2, exp_factor);
    lapl_basis_functions[20] = new lapl_HarmonicOscillator_20(k, k2, exp_factor);
    lapl_basis_functions[21] = new lapl_HarmonicOscillator_21(k, k2, exp_factor);
    lapl_basis_functions[22] = new lapl_HarmonicOscillator_22(k, k2, exp_factor);
    lapl_basis_functions[23] = new lapl_HarmonicOscillator_23(k, k2, exp_factor);
    lapl_basis_functions[24] = new lapl_HarmonicOscillator_24(k, k2, exp_factor);
    lapl_basis_functions[25] = new lapl_HarmonicOscillator_25(k, k2, exp_factor);
    lapl_basis_functions[26] = new lapl_HarmonicOscillator_26(k, k2, exp_factor);
    lapl_basis_functions[27] = new lapl_HarmonicOscillator_27(k, k2, exp_factor);

}

void AlphaHarmonicOscillator::get_qnums() {
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
    } else if (n == 5) {
        double x2 = x*x;
        return 32 * x2 * x2 * x - 160 * x2 * x + 120 * x;
    } else if (n == 6) {
        double x2 = x*x;
        double x4 = x2*x2;
        return 64 * x4 * x2 - 480 * x4 + 720 * x2 - 120;
    } else {
        std::cout << "Unsopported Hermite polynomial level: " << n << std::endl;
        exit(1);
    }
}

double AlphaHarmonicOscillator::get_dell_alpha_phi(Walker* walker, int p, int q_num) {

    int nij;
    double rij;

    double H_fac = 0;
    for (int j = 0; j < dim; j++) {

        rij = walker->r(p, j);
        nij = qnums(q_num, j);

        H_fac += rij * nij * H(nij - 1, (*k) * rij) / H(nij, (*k) * rij);
    }

    H_fac *= (*k) / (*alpha);

    return H_fac - 0.5 * w * walker->get_r_i2(p);

}

void AlphaHarmonicOscillator::set_qnum_indie_terms(Walker* walker, int i) {
    *exp_factor = exp(-0.5 * (*k2) * walker->get_r_i2(i));
}

double AlphaHarmonicOscillator::get_coulomb_element(const arma::uvec& qnum_set) {
    //Needs transformations to be effective. Simen Kvaal style.
    //    using namespace arma;
    //
    //    int n_x, n_y;
    //
    //    ivec n_set(4);
    //    ivec m_set(4);
    //
    //    for (int i = 0; i < 4; i++) {
    //        n_x = qnums(qnum_set(i), 0);
    //        n_y = qnums(qnum_set(i), 1);
    //
    //        m_set(i) = n_x - n_y;
    //        n_set(i) = (n_x + n_y - abs(m_set(i))) / 2;
    //
    //    }
    //    //    cout << "q " <<qnum_set.st() << "m " << m_set.st() << "n " << n_set.st() << endl;
    ////    double element = coulomb(n_set(0), m_set(0),
    ////            n_set(1), m_set(1),
    ////            n_set(3), m_set(3),
    ////            n_set(2), m_set(2));
    //    //    cout << element << "\n-------" << endl;
    //
    //    int n_c = 10000000;
    //    double a = -3;
    //    double b = 3;
    //    Diffusion* diff = new Simple(1, 1, 1, 1000, 1);
    //
    //    Walker* dummy = new Walker(2, dim);
    //    double I = 0;
    //    for (int i = 0; i < n_c; i++) {
    //        dummy->r = a + (b - a) * randu(2, dim);
    //
    //        set_qnum_indie_terms(dummy, 0);
    //        double local_I = phi(dummy, 0, qnum_set(0)) * phi(dummy, 0, qnum_set(2));
    //
    //        set_qnum_indie_terms(dummy, 1);
    //        local_I *= phi(dummy, 1, qnum_set(1)) * phi(dummy, 1, qnum_set(3));
    //
    //        I += local_I / dummy->calc_r_rel(0, 1);
    //
    //
    //    }
    //
    //    I *= (b - a) / n_c;
    //    if (I < 1E-2) {
    //        I = 0;
    //    }
    //    //    cout << I << "  " << element << endl;
    //
    //    return I;

    (void) qnum_set;
    return 0;

}

double AlphaHarmonicOscillator::get_sp_energy(int qnum) const {
    int n_x = qnums(qnum, 0);
    int n_y = qnums(qnum, 1);

    return w * (n_x + n_y + 1);

}
