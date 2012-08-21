/* 
 * File:   Orbitals.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

#include "../QMCheaders.h"

Orbitals::Orbitals(int n_p, int dim) {
    this->n_p = n_p;
    this->dim = dim;

    int max_implemented = 15; //for 30 particles
    basis_functions = new function*[max_implemented];
}

oscillator_basis::oscillator_basis(int n_p, int dim, double alpha, double w)
: Orbitals(n_p, dim) {
    this->alpha = alpha;
    this->w = w;

    int max_implemented = 15; //for 30 particles

    dell_basis_functions = new function**[dim];
    for (int i = 0; i < dim; i++) {
        dell_basis_functions[i] = new function*[max_implemented];
    }

    lapl_basis_functions = new function*[max_implemented];


    basis_functions[0] = new HO_1(alpha, w);
    basis_functions[1] = new HO_2(alpha, w);
    basis_functions[2] = new HO_3(alpha, w);
    basis_functions[3] = new HO_4(alpha, w);
    basis_functions[4] = new HO_5(alpha, w);
    basis_functions[5] = new HO_6(alpha, w);
    basis_functions[6] = new HO_7(alpha, w);
    basis_functions[7] = new HO_8(alpha, w);
    basis_functions[8] = new HO_9(alpha, w);
    basis_functions[9] = new HO_10(alpha, w);
    basis_functions[10] = new HO_11(alpha, w);
    basis_functions[11] = new HO_12(alpha, w);
    basis_functions[12] = new HO_13(alpha, w);
    basis_functions[13] = new HO_14(alpha, w);
    basis_functions[14] = new HO_15(alpha, w);

    dell_basis_functions[0][0] = new dell_HO_1_0(alpha, w);
    dell_basis_functions[1][0] = new dell_HO_1_1(alpha, w);
    dell_basis_functions[0][1] = new dell_HO_2_0(alpha, w);
    dell_basis_functions[1][1] = new dell_HO_2_1(alpha, w);
    dell_basis_functions[0][2] = new dell_HO_3_0(alpha, w);
    dell_basis_functions[1][2] = new dell_HO_3_1(alpha, w);
    dell_basis_functions[0][3] = new dell_HO_4_0(alpha, w);
    dell_basis_functions[1][3] = new dell_HO_4_1(alpha, w);
    dell_basis_functions[0][4] = new dell_HO_5_0(alpha, w);
    dell_basis_functions[1][4] = new dell_HO_5_1(alpha, w);
    dell_basis_functions[0][5] = new dell_HO_6_0(alpha, w);
    dell_basis_functions[1][5] = new dell_HO_6_1(alpha, w);
    dell_basis_functions[0][6] = new dell_HO_7_0(alpha, w);
    dell_basis_functions[1][6] = new dell_HO_7_1(alpha, w);
    dell_basis_functions[0][7] = new dell_HO_8_0(alpha, w);
    dell_basis_functions[1][7] = new dell_HO_8_1(alpha, w);
    dell_basis_functions[0][8] = new dell_HO_9_0(alpha, w);
    dell_basis_functions[1][8] = new dell_HO_9_1(alpha, w);
    dell_basis_functions[0][9] = new dell_HO_10_0(alpha, w);
    dell_basis_functions[1][9] = new dell_HO_10_1(alpha, w);
    dell_basis_functions[0][10] = new dell_HO_11_0(alpha, w);
    dell_basis_functions[1][10] = new dell_HO_11_1(alpha, w);
    dell_basis_functions[0][11] = new dell_HO_12_0(alpha, w);
    dell_basis_functions[1][11] = new dell_HO_12_1(alpha, w);
    dell_basis_functions[0][12] = new dell_HO_13_0(alpha, w);
    dell_basis_functions[1][12] = new dell_HO_13_1(alpha, w);
    dell_basis_functions[0][13] = new dell_HO_14_0(alpha, w);
    dell_basis_functions[1][13] = new dell_HO_14_1(alpha, w);
    dell_basis_functions[0][14] = new dell_HO_15_0(alpha, w);
    dell_basis_functions[1][14] = new dell_HO_15_1(alpha, w);


    lapl_basis_functions[0] = new lapl_HO_1(alpha, w);
    lapl_basis_functions[1] = new lapl_HO_2(alpha, w);
    lapl_basis_functions[2] = new lapl_HO_3(alpha, w);
    lapl_basis_functions[3] = new lapl_HO_4(alpha, w);
    lapl_basis_functions[4] = new lapl_HO_5(alpha, w);
    lapl_basis_functions[5] = new lapl_HO_6(alpha, w);
    lapl_basis_functions[6] = new lapl_HO_7(alpha, w);
    lapl_basis_functions[7] = new lapl_HO_8(alpha, w);
    lapl_basis_functions[8] = new lapl_HO_9(alpha, w);
    lapl_basis_functions[9] = new lapl_HO_10(alpha, w);
    lapl_basis_functions[10] = new lapl_HO_11(alpha, w);
    lapl_basis_functions[11] = new lapl_HO_12(alpha, w);
    lapl_basis_functions[12] = new lapl_HO_13(alpha, w);
    lapl_basis_functions[13] = new lapl_HO_14(alpha, w);
    lapl_basis_functions[14] = new lapl_HO_15(alpha, w);




}

double oscillator_basis::phi(const Walker* walker, int particle, int q_num) const {
    return basis_functions[q_num]->eval(walker, particle);
}

double oscillator_basis::del_phi(const Walker* walker, int particle, int q_num, int d) const {
    return dell_basis_functions[d][q_num]->eval(walker, particle);
}

double oscillator_basis::lapl_phi(const Walker* walker, int particle, int q_num) const {
    return lapl_basis_functions[q_num]->eval(walker, particle);
}

