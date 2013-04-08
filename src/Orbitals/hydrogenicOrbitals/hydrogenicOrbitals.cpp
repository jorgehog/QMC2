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

}

double hydrogenicOrbitals::get_dell_alpha_phi(const Walker* walker, int qnum, int i) {

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

double hydrogenicOrbitals::get_coulomb_element(const arma::uvec& qnum_set) {

    int Mleft = qnums(qnum_set(0), 2) + qnums(qnum_set(1), 2);
    int Mright = qnums(qnum_set(2), 2) + qnums(qnum_set(3), 2);

    if (Mleft != Mright) {
        return 0;
    }
//
//    arma::vec N(4);
//    N.zeros();
//
//    int n_c = 1000000;
//    int m = qnum_set.min();
//    //    double L = -0.2*(Z<11) - 0.2*(Z < 5) + 0.6*(Z < 3);
//    double L = 5;
//    double a = -L;
//    double b = L;
//
//    Walker* dummy = new Walker(2, dim);
//    double I = 0;
//    double r_rel;
//    int x = 0;
//
//    for (int i = 0; i < n_c; i++) {
//
//        dummy->r = a + (b - a) * arma::randu(2, dim);
//        r_rel = dummy->calc_r_rel(0, 1);
//        if (r_rel < 1E-3) {
//            x++;
//            continue;
//        }
//        dummy->calc_r_i2();
//        set_qnum_indie_terms(dummy, 0);
//        double local_I = phi(dummy, 0, qnum_set(0)) * phi(dummy, 0, qnum_set(2));
//
//        set_qnum_indie_terms(dummy, 1);
//        local_I *= phi(dummy, 1, qnum_set(1)) * phi(dummy, 1, qnum_set(3));
//
//        I += local_I / r_rel;
//
//
//    }
//
//    I *= pow((b - a), 2 * dim) / (n_c - x);
//    //    double L = -0.5*(Z<11) - 1*(Z < 5) + 2*(Z < 3);
//
//    L = 3;
//    a = -L;
//    b = L;
//
//    for (int i = 0; i < n_c; i++) {
//
//        dummy->r = a + (b - a) * arma::randu(2, dim);
//        dummy->calc_r_i2();
//        set_qnum_indie_terms(dummy, 0);
//
//        N(0) += phi(dummy, 0, qnum_set(0)) * phi(dummy, 0, qnum_set(0));
//        N(2) += phi(dummy, 0, qnum_set(2)) * phi(dummy, 0, qnum_set(2));
//
//        set_qnum_indie_terms(dummy, 1);
//
//        N(1) += phi(dummy, 1, qnum_set(1)) * phi(dummy, 1, qnum_set(1));
//        N(3) += phi(dummy, 1, qnum_set(3)) * phi(dummy, 1, qnum_set(3));
//
//    }
//
//    N *= pow((b - a), dim) / n_c;
//    N = arma::ones<arma::vec > (4) / arma::sqrt(N);
//    std::cout << N << std::endl; a = -L;
//    b = L;
//
//    for (int i = 0; i < n_c; i++) {
//
//        dummy->r = a + (b - a) * arma::randu(2, dim);
//        dummy->calc_r_i2();
//        set_qnum_indie_terms(dummy, 0);
//
//        N(0) += phi(dummy, 0, qnum_set(0)) * phi(dummy, 0, qnum_set(0));
//        N(2) += phi(dummy, 0, qnum_set(2)) * phi(dummy, 0, qnum_set(2));
//
//        set_qnum_indie_terms(dummy, 1);
//
//        N(1) += phi(dummy, 1, qnum_set(1)) * phi(dummy, 1, qnum_set(1));
//        N(3) += phi(dummy, 1, qnum_set(3)) * phi(dummy, 1, qnum_set(3));
//
//    }
//
//    N *= pow((b - a), dim) / n_c;
//    N = arma::ones<arma::vec > (4) / arma::sqrt(N);
//    std::cout << N << std::endl;
//

    //    std::cout << qnum_set << std::endl;
    //    std::cout << I << std::endl;
    //    return 0;

    int n1 = qnums(qnum_set(0), 0);
    int n2 = qnums(qnum_set(1), 0);
    int n3 = qnums(qnum_set(2), 0);
    int n4 = qnums(qnum_set(3), 0);

    using namespace arma;

    if (n1 == 1 && n2 == 1 && n3 == 1 && n4 == 1) return (5 * Z) / 8.0;
    else if (n1 == 1 && n2 == 1 && n3 == 1 && n4 == 2) return (4096 * sqrt(2) * Z) / 64827.0;
    else if (n1 == 1 && n2 == 1 && n3 == 1 && n4 == 3) return (1269 * sqrt(3) * Z) / 50000.0;
    else if (n1 == 1 && n2 == 1 && n3 == 1 && n4 == 4) return (416415744 * Z) / 15083778125.0;
    else if (n1 == 1 && n2 == 1 && n3 == 2 && n4 == 1) return (4096 * sqrt(2) * Z) / 64827.0;
    else if (n1 == 1 && n2 == 1 && n3 == 2 && n4 == 2) return (16 * Z) / 729.0;
    else if (n1 == 1 && n2 == 1 && n3 == 2 && n4 == 3) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 1 && n2 == 1 && n3 == 2 && n4 == 4) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 1 && n2 == 1 && n3 == 3 && n4 == 1) return (1269 * sqrt(3) * Z) / 50000.0;
    else if (n1 == 1 && n2 == 1 && n3 == 3 && n4 == 2) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 1 && n2 == 1 && n3 == 3 && n4 == 3) return (189 * Z) / 32768.0;
    else if (n1 == 1 && n2 == 1 && n3 == 3 && n4 == 4) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 1 && n2 == 1 && n3 == 4 && n4 == 1) return (416415744 * Z) / 15083778125.0;
    else if (n1 == 1 && n2 == 1 && n3 == 4 && n4 == 2) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 1 && n2 == 1 && n3 == 4 && n4 == 3) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 1 && n2 == 1 && n3 == 4 && n4 == 4) return (22848 * Z) / 9765625.0;
    else if (n1 == 1 && n2 == 2 && n3 == 1 && n4 == 1) return (4096 * sqrt(2) * Z) / 64827.0;
    else if (n1 == 1 && n2 == 2 && n3 == 1 && n4 == 2) return (17 * Z) / 81.0;
    else if (n1 == 1 && n2 == 2 && n3 == 1 && n4 == 3) return (1555918848 * sqrt(6) * Z) / 75429903125.0;
    else if (n1 == 1 && n2 == 2 && n3 == 1 && n4 == 4) return (288456704 * sqrt(2) * Z) / 14206147659.0;
    else if (n1 == 1 && n2 == 2 && n3 == 2 && n4 == 1) return (16 * Z) / 729.0;
    else if (n1 == 1 && n2 == 2 && n3 == 2 && n4 == 2) return (512 * sqrt(2) * Z) / 84375.0;
    else if (n1 == 1 && n2 == 2 && n3 == 2 && n4 == 3) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 1 && n2 == 2 && n3 == 2 && n4 == 4) return (376832 * Z) / 129140163.0;
    else if (n1 == 1 && n2 == 2 && n3 == 3 && n4 == 1) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 1 && n2 == 2 && n3 == 3 && n4 == 2) return (29943 * sqrt(3) * Z) / 13176688.0;
    else if (n1 == 1 && n2 == 2 && n3 == 3 && n4 == 3) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 1 && n2 == 2 && n3 == 3 && n4 == 4) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 1 && n2 == 2 && n3 == 4 && n4 == 1) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 1 && n2 == 2 && n3 == 4 && n4 == 2) return (108838912 * Z) / 44840334375.0;
    else if (n1 == 1 && n2 == 2 && n3 == 4 && n4 == 3) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 1 && n2 == 2 && n3 == 4 && n4 == 4) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 1 && n2 == 3 && n3 == 1 && n4 == 1) return (1269 * sqrt(3) * Z) / 50000.0;
    else if (n1 == 1 && n2 == 3 && n3 == 1 && n4 == 2) return (1555918848 * sqrt(6) * Z) / 75429903125.0;
    else if (n1 == 1 && n2 == 3 && n3 == 1 && n4 == 3) return (815 * Z) / 8192.0;
    else if (n1 == 1 && n2 == 3 && n3 == 1 && n4 == 4) return (11694850770862080 * sqrt(3) * Z) / 702392443647273463.0;
    else if (n1 == 1 && n2 == 3 && n3 == 2 && n4 == 1) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 1 && n2 == 3 && n3 == 2 && n4 == 2) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 1 && n2 == 3 && n3 == 2 && n4 == 3) return (37826560 * sqrt(2) * Z) / 22024729467.0;
    else if (n1 == 1 && n2 == 3 && n3 == 2 && n4 == 4) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 1 && n2 == 3 && n3 == 3 && n4 == 1) return (189 * Z) / 32768.0;
    else if (n1 == 1 && n2 == 3 && n3 == 3 && n4 == 2) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 1 && n2 == 3 && n3 == 3 && n4 == 3) return (617 * Z) / (314928 * sqrt(3));
    else if (n1 == 1 && n2 == 3 && n3 == 3 && n4 == 4) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 1 && n2 == 3 && n3 == 4 && n4 == 1) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 1 && n2 == 3 && n3 == 4 && n4 == 2) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 1 && n2 == 3 && n3 == 4 && n4 == 3) return (90581886173184 * Z) / 129457847542653125.0;
    else if (n1 == 1 && n2 == 3 && n3 == 4 && n4 == 4) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 1 && n2 == 4 && n3 == 1 && n4 == 1) return (416415744 * Z) / 15083778125.0;
    else if (n1 == 1 && n2 == 4 && n3 == 1 && n4 == 2) return (288456704 * sqrt(2) * Z) / 14206147659.0;
    else if (n1 == 1 && n2 == 4 && n3 == 1 && n4 == 3) return (11694850770862080 * sqrt(3) * Z) / 702392443647273463.0;
    else if (n1 == 1 && n2 == 4 && n3 == 1 && n4 == 4) return (22513 * Z) / 390625.0;
    else if (n1 == 1 && n2 == 4 && n3 == 2 && n4 == 1) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 1 && n2 == 4 && n3 == 2 && n4 == 2) return (376832 * Z) / 129140163.0;
    else if (n1 == 1 && n2 == 4 && n3 == 2 && n4 == 3) return (487489536 * sqrt(6) * Z) / 7629394531255.0;
    else if (n1 == 1 && n2 == 4 && n3 == 2 && n4 == 4) return (5053 * Z) / (3538944 * sqrt(2));
    else if (n1 == 1 && n2 == 4 && n3 == 3 && n4 == 1) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 1 && n2 == 4 && n3 == 3 && n4 == 2) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 1 && n2 == 4 && n3 == 3 && n4 == 3) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 1 && n2 == 4 && n3 == 3 && n4 == 4) return (1243165779 * sqrt(3) * Z) / 4564986729776.0;
    else if (n1 == 1 && n2 == 4 && n3 == 4 && n4 == 1) return (22848 * Z) / 9765625.0;
    else if (n1 == 1 && n2 == 4 && n3 == 4 && n4 == 2) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 1 && n2 == 4 && n3 == 4 && n4 == 3) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 1 && n2 == 4 && n3 == 4 && n4 == 4) return (1804351488 * Z) / 6179146071875.0;
    else if (n1 == 2 && n2 == 1 && n3 == 1 && n4 == 1) return (4096 * sqrt(2) * Z) / 64827.0;
    else if (n1 == 2 && n2 == 1 && n3 == 1 && n4 == 2) return (16 * Z) / 729.0;
    else if (n1 == 2 && n2 == 1 && n3 == 1 && n4 == 3) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 2 && n2 == 1 && n3 == 1 && n4 == 4) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 2 && n2 == 1 && n3 == 2 && n4 == 1) return (17 * Z) / 81.0;
    else if (n1 == 2 && n2 == 1 && n3 == 2 && n4 == 2) return (512 * sqrt(2) * Z) / 84375.0;
    else if (n1 == 2 && n2 == 1 && n3 == 2 && n4 == 3) return (29943 * sqrt(3) * Z) / 13176688.0;
    else if (n1 == 2 && n2 == 1 && n3 == 2 && n4 == 4) return (108838912 * Z) / 44840334375.0;
    else if (n1 == 2 && n2 == 1 && n3 == 3 && n4 == 1) return (1555918848 * sqrt(6) * Z) / 75429903125.0;
    else if (n1 == 2 && n2 == 1 && n3 == 3 && n4 == 2) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 2 && n2 == 1 && n3 == 3 && n4 == 3) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 2 && n2 == 1 && n3 == 3 && n4 == 4) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 2 && n2 == 1 && n3 == 4 && n4 == 1) return (288456704 * sqrt(2) * Z) / 14206147659.0;
    else if (n1 == 2 && n2 == 1 && n3 == 4 && n4 == 2) return (376832 * Z) / 129140163.0;
    else if (n1 == 2 && n2 == 1 && n3 == 4 && n4 == 3) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 2 && n2 == 1 && n3 == 4 && n4 == 4) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 2 && n2 == 2 && n3 == 1 && n4 == 1) return (16 * Z) / 729.0;
    else if (n1 == 2 && n2 == 2 && n3 == 1 && n4 == 2) return (512 * sqrt(2) * Z) / 84375.0;
    else if (n1 == 2 && n2 == 2 && n3 == 1 && n4 == 3) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 2 && n2 == 2 && n3 == 1 && n4 == 4) return (376832 * Z) / 129140163.0;
    else if (n1 == 2 && n2 == 2 && n3 == 2 && n4 == 1) return (512 * sqrt(2) * Z) / 84375.0;
    else if (n1 == 2 && n2 == 2 && n3 == 2 && n4 == 2) return (77 * Z) / 512.0;
    else if (n1 == 2 && n2 == 2 && n3 == 2 && n4 == 3) return (5870679552 * sqrt(6) * Z) / 669871503125.0;
    else if (n1 == 2 && n2 == 2 && n3 == 2 && n4 == 4) return (31363072 * sqrt(2) * Z) / 4202539929.0;
    else if (n1 == 2 && n2 == 2 && n3 == 3 && n4 == 1) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 2 && n2 == 2 && n3 == 3 && n4 == 2) return (5870679552 * sqrt(6) * Z) / 669871503125.0;
    else if (n1 == 2 && n2 == 2 && n3 == 3 && n4 == 3) return (73008 * Z) / 9765625.0;
    else if (n1 == 2 && n2 == 2 && n3 == 3 && n4 == 4) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 2 && n2 == 2 && n3 == 4 && n4 == 1) return (376832 * Z) / 129140163.0;
    else if (n1 == 2 && n2 == 2 && n3 == 4 && n4 == 2) return (31363072 * sqrt(2) * Z) / 4202539929.0;
    else if (n1 == 2 && n2 == 2 && n3 == 4 && n4 == 3) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 2 && n2 == 2 && n3 == 4 && n4 == 4) return (424 * Z) / 177147.0;
    else if (n1 == 2 && n2 == 3 && n3 == 1 && n4 == 1) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 2 && n2 == 3 && n3 == 1 && n4 == 2) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 2 && n2 == 3 && n3 == 1 && n4 == 3) return (37826560 * sqrt(2) * Z) / 22024729467.0;
    else if (n1 == 2 && n2 == 3 && n3 == 1 && n4 == 4) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 2 && n2 == 3 && n3 == 2 && n4 == 1) return (29943 * sqrt(3) * Z) / 13176688.0;
    else if (n1 == 2 && n2 == 3 && n3 == 2 && n4 == 2) return (5870679552 * sqrt(6) * Z) / 669871503125.0;
    else if (n1 == 2 && n2 == 3 && n3 == 2 && n4 == 3) return (32857 * Z) / 390625.0;
    else if (n1 == 2 && n2 == 3 && n3 == 2 && n4 == 4) return (55508689880137728 * sqrt(3) * Z) / 5049196699148208943.0;
    else if (n1 == 2 && n2 == 3 && n3 == 3 && n4 == 1) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 2 && n2 == 3 && n3 == 3 && n4 == 2) return (73008 * Z) / 9765625.0;
    else if (n1 == 2 && n2 == 3 && n3 == 3 && n4 == 3) return (6890942464 * sqrt(2 / 3.0) * Z) / 1210689028125.0;
    else if (n1 == 2 && n2 == 3 && n3 == 3 && n4 == 4) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 2 && n2 == 3 && n3 == 4 && n4 == 1) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 2 && n2 == 3 && n3 == 4 && n4 == 2) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 2 && n2 == 3 && n3 == 4 && n4 == 3) return (36645380390912 * sqrt(2) * Z) / 24984212408264457.0;
    else if (n1 == 2 && n2 == 3 && n3 == 4 && n4 == 4) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 2 && n2 == 4 && n3 == 1 && n4 == 1) return (98304 * sqrt(2) * Z) / 194871716.0;
    else if (n1 == 2 && n2 == 4 && n3 == 1 && n4 == 2) return (376832 * Z) / 129140163.0;
    else if (n1 == 2 && n2 == 4 && n3 == 1 && n4 == 3) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 2 && n2 == 4 && n3 == 1 && n4 == 4) return (5053 * Z) / (3538944 * sqrt(2));
    else if (n1 == 2 && n2 == 4 && n3 == 2 && n4 == 1) return (108838912 * Z) / 44840334375.0;
    else if (n1 == 2 && n2 == 4 && n3 == 2 && n4 == 2) return (31363072 * sqrt(2) * Z) / 4202539929.0;
    else if (n1 == 2 && n2 == 4 && n3 == 2 && n4 == 3) return (55508689880137728 * sqrt(3) * Z) / 5049196699148208943.0;
    else if (n1 == 2 && n2 == 4 && n3 == 2 && n4 == 4) return (4043 * Z) / 78732.0;
    else if (n1 == 2 && n2 == 4 && n3 == 3 && n4 == 1) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 2 && n2 == 4 && n3 == 3 && n4 == 2) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 2 && n2 == 4 && n3 == 3 && n4 == 3) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 2 && n2 == 4 && n3 == 3 && n4 == 4) return (2496169683 * sqrt(1.5) * Z) / 1677721600000.0;
    else if (n1 == 2 && n2 == 4 && n3 == 4 && n4 == 1) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 2 && n2 == 4 && n3 == 4 && n4 == 2) return (424 * Z) / 177147.0;
    else if (n1 == 2 && n2 == 4 && n3 == 4 && n4 == 3) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 2 && n2 == 4 && n3 == 4 && n4 == 4) return (21252608 * sqrt(2) * Z) / 35595703125.0;
    else if (n1 == 3 && n2 == 1 && n3 == 1 && n4 == 1) return (1269 * sqrt(3) * Z) / 50000.0;
    else if (n1 == 3 && n2 == 1 && n3 == 1 && n4 == 2) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 3 && n2 == 1 && n3 == 1 && n4 == 3) return (189 * Z) / 32768.0;
    else if (n1 == 3 && n2 == 1 && n3 == 1 && n4 == 4) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 3 && n2 == 1 && n3 == 2 && n4 == 1) return (1555918848 * sqrt(6) * Z) / 75429903125.0;
    else if (n1 == 3 && n2 == 1 && n3 == 2 && n4 == 2) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 3 && n2 == 1 && n3 == 2 && n4 == 3) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 3 && n2 == 1 && n3 == 2 && n4 == 4) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 3 && n2 == 1 && n3 == 3 && n4 == 1) return (815 * Z) / 8192.0;
    else if (n1 == 3 && n2 == 1 && n3 == 3 && n4 == 2) return (37826560 * sqrt(2) * Z) / 22024729467.0;
    else if (n1 == 3 && n2 == 1 && n3 == 3 && n4 == 3) return (617 * Z) / (314928 * sqrt(3));
    else if (n1 == 3 && n2 == 1 && n3 == 3 && n4 == 4) return (90581886173184 * Z) / 129457847542653125.0;
    else if (n1 == 3 && n2 == 1 && n3 == 4 && n4 == 1) return (11694850770862080 * sqrt(3) * Z) / 702392443647273463.0;
    else if (n1 == 3 && n2 == 1 && n3 == 4 && n4 == 2) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 3 && n2 == 1 && n3 == 4 && n4 == 3) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 3 && n2 == 1 && n3 == 4 && n4 == 4) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 3 && n2 == 2 && n3 == 1 && n4 == 1) return (110592 * sqrt(6) * Z) / 24137569.0;
    else if (n1 == 3 && n2 == 2 && n3 == 1 && n4 == 2) return (29943 * sqrt(3) * Z) / 13176688.0;
    else if (n1 == 3 && n2 == 2 && n3 == 1 && n4 == 3) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 3 && n2 == 2 && n3 == 1 && n4 == 4) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 3 && n2 == 2 && n3 == 2 && n4 == 1) return (2160 * sqrt(3) * Z) / 823543.0;
    else if (n1 == 3 && n2 == 2 && n3 == 2 && n4 == 2) return (5870679552 * sqrt(6) * Z) / 669871503125.0;
    else if (n1 == 3 && n2 == 2 && n3 == 2 && n4 == 3) return (73008 * Z) / 9765625.0;
    else if (n1 == 3 && n2 == 2 && n3 == 2 && n4 == 4) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 3 && n2 == 2 && n3 == 3 && n4 == 1) return (37826560 * sqrt(2) * Z) / 22024729467.0;
    else if (n1 == 3 && n2 == 2 && n3 == 3 && n4 == 2) return (32857 * Z) / 390625.0;
    else if (n1 == 3 && n2 == 2 && n3 == 3 && n4 == 3) return (6890942464 * sqrt(2 / 3.0) * Z) / 1210689028125.0;
    else if (n1 == 3 && n2 == 2 && n3 == 3 && n4 == 4) return (36645380390912 * sqrt(2) * Z) / 24984212408264457.0;
    else if (n1 == 3 && n2 == 2 && n3 == 4 && n4 == 1) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 3 && n2 == 2 && n3 == 4 && n4 == 2) return (55508689880137728 * sqrt(3) * Z) / 5049196699148208943.0;
    else if (n1 == 3 && n2 == 2 && n3 == 4 && n4 == 3) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 3 && n2 == 2 && n3 == 4 && n4 == 4) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 3 && n2 == 3 && n3 == 1 && n4 == 1) return (189 * Z) / 32768.0;
    else if (n1 == 3 && n2 == 3 && n3 == 1 && n4 == 2) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 3 && n2 == 3 && n3 == 1 && n4 == 3) return (617 * Z) / (314928 * sqrt(3));
    else if (n1 == 3 && n2 == 3 && n3 == 1 && n4 == 4) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 3 && n2 == 3 && n3 == 2 && n4 == 1) return (1216512 * sqrt(2) * Z) / 815730721.0;
    else if (n1 == 3 && n2 == 3 && n3 == 2 && n4 == 2) return (73008 * Z) / 9765625.0;
    else if (n1 == 3 && n2 == 3 && n3 == 2 && n4 == 3) return (6890942464 * sqrt(2/3.0) * Z) / 1210689028125.0;
    else if (n1 == 3 && n2 == 3 && n3 == 2 && n4 == 4) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 3 && n2 == 3 && n3 == 3 && n4 == 1) return (617 * Z) / (314928 * sqrt(3));
    else if (n1 == 3 && n2 == 3 && n3 == 3 && n4 == 2) return (6890942464 * sqrt(2/3.0) * Z) / 1210689028125.0;
    else if (n1 == 3 && n2 == 3 && n3 == 3 && n4 == 3) return (17 * Z) / 2567.0;
    else if (n1 == 3 && n2 == 3 && n3 == 3 && n4 == 4) return (2486755845603328 * Z) / (158298797548828125 * sqrt(3));
    else if (n1 == 3 && n2 == 3 && n3 == 4 && n4 == 1) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 3 && n2 == 3 && n3 == 4 && n4 == 2) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 3 && n2 == 3 && n3 == 4 && n4 == 3) return (2486755845603328 * Z) / (158298797548828125 * sqrt(3));
    else if (n1 == 3 && n2 == 3 && n3 == 4 && n4 == 4) return (2560158144 * Z) / 678223072849.0;
    else if (n1 == 3 && n2 == 4 && n3 == 1 && n4 == 1) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 3 && n2 == 4 && n3 == 1 && n4 == 2) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 3 && n2 == 4 && n3 == 1 && n4 == 3) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 3 && n2 == 4 && n3 == 1 && n4 == 4) return (1243165779 * sqrt(3) * Z) / 4564986729776.0;
    else if (n1 == 3 && n2 == 4 && n3 == 2 && n4 == 1) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 3 && n2 == 4 && n3 == 2 && n4 == 2) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 3 && n2 == 4 && n3 == 2 && n4 == 3) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 3 && n2 == 4 && n3 == 2 && n4 == 4) return (2496169683 * sqrt(1.5) * Z) / 1677721600000.0;
    else if (n1 == 3 && n2 == 4 && n3 == 3 && n4 == 1) return (90581886173184 * Z) / 129457847542653125.0;
    else if (n1 == 3 && n2 == 4 && n3 == 3 && n4 == 2) return (36645380390912 * sqrt(2) * Z) / 24984212408264457.0;
    else if (n1 == 3 && n2 == 4 && n3 == 3 && n4 == 3) return (2486755845603328 * Z) / (158298797548828125 * sqrt(3));
    else if (n1 == 3 && n2 == 4 && n3 == 3 && n4 == 4) return (621550729 * Z) / 13841287201.0;
    else if (n1 == 3 && n2 == 4 && n3 == 4 && n4 == 1) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 3 && n2 == 4 && n3 == 4 && n4 == 2) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 3 && n2 == 4 && n3 == 4 && n4 == 3) return (2560158144 * Z) / 678223072849.0;
    else if (n1 == 3 && n2 == 4 && n3 == 4 && n4 == 4) return (413631006610176000 * sqrt(3) * Z) / 249430673908303812379.0;
    else if (n1 == 4 && n2 == 1 && n3 == 1 && n4 == 1) return (416415744 * Z) / 15083778125.0;
    else if (n1 == 4 && n2 == 1 && n3 == 1 && n4 == 2) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 4 && n2 == 1 && n3 == 1 && n4 == 3) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 4 && n2 == 1 && n3 == 1 && n4 == 4) return (22848 * Z) / 9765625.0;
    else if (n1 == 4 && n2 == 1 && n3 == 2 && n4 == 1) return (288456704 * sqrt(2) * Z) / 14206147659.0;
    else if (n1 == 4 && n2 == 1 && n3 == 2 && n4 == 2) return (376832 * Z) / 129140163.0;
    else if (n1 == 4 && n2 == 1 && n3 == 2 && n4 == 3) return (423788544 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 4 && n2 == 1 && n3 == 2 && n4 == 4) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 4 && n2 == 1 && n3 == 3 && n4 == 1) return (11694850770862080 * sqrt(3) * Z) / 702392443647273463.0;
    else if (n1 == 4 && n2 == 1 && n3 == 3 && n4 == 2) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 4 && n2 == 1 && n3 == 3 && n4 == 3) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 4 && n2 == 1 && n3 == 3 && n4 == 4) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 4 && n2 == 1 && n3 == 4 && n4 == 1) return (22513 * Z) / 390625.0;
    else if (n1 == 4 && n2 == 1 && n3 == 4 && n4 == 2) return (5053 * Z) / (3538944 * sqrt(2));
    else if (n1 == 4 && n2 == 1 && n3 == 4 && n4 == 3) return (1243165779 * sqrt(3) * Z) / 4564986729776.0;
    else if (n1 == 4 && n2 == 1 && n3 == 4 && n4 == 4) return (1804351488 * Z) / 6179146071875.0;
    else if (n1 == 4 && n2 == 2 && n3 == 1 && n4 == 1) return (98304 * sqrt(2) * Z) / 19487171.0;
    else if (n1 == 4 && n2 == 2 && n3 == 1 && n4 == 2) return (108838912 * Z) / 44840334375.0;
    else if (n1 == 4 && n2 == 2 && n3 == 1 && n4 == 3) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 4 && n2 == 2 && n3 == 1 && n4 == 4) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 4 && n2 == 2 && n3 == 2 && n4 == 1) return (376832 * Z) / 129140163.0;
    else if (n1 == 4 && n2 == 2 && n3 == 2 && n4 == 2) return (31363072 * sqrt(2) * Z) / 4202539929.0;
    else if (n1 == 4 && n2 == 2 && n3 == 2 && n4 == 3) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 4 && n2 == 2 && n3 == 2 && n4 == 4) return (424 * Z) / 177147.0;
    else if (n1 == 4 && n2 == 2 && n3 == 3 && n4 == 1) return (487489536 * sqrt(6) * Z) / 762939453125.0;
    else if (n1 == 4 && n2 == 2 && n3 == 3 && n4 == 2) return (55508689880137728 * sqrt(3) * Z) / 5049196699148208943.0;
    else if (n1 == 4 && n2 == 2 && n3 == 3 && n4 == 3) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 4 && n2 == 2 && n3 == 3 && n4 == 4) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 4 && n2 == 2 && n3 == 4 && n4 == 1) return (5053 * Z) / (3538944 * sqrt(2));
    else if (n1 == 4 && n2 == 2 && n3 == 4 && n4 == 2) return (4043 * Z) / 78732.0;
    else if (n1 == 4 && n2 == 2 && n3 == 4 && n4 == 3) return (2496169683 * sqrt(1.5) * Z) / 1677721600000.0;
    else if (n1 == 4 && n2 == 2 && n3 == 4 && n4 == 4) return (21252608 * sqrt(2) * Z) / 35595703125.0;
    else if (n1 == 4 && n2 == 3 && n3 == 1 && n4 == 1) return (1808400384 * sqrt(3) * Z) / 852891037441.0;
    else if (n1 == 4 && n2 == 3 && n3 == 1 && n4 == 2) return (10148806656 * sqrt(6) * Z) / 19073486328125.0;
    else if (n1 == 4 && n2 == 3 && n3 == 1 && n4 == 3) return (90581886173184 * Z) / 129457847542653125.0;
    else if (n1 == 4 && n2 == 3 && n3 == 1 && n4 == 4) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 4 && n2 == 3 && n3 == 2 && n4 == 1) return (423788544 * sqrt(6) * Z) / 7629394531258.0;
    else if (n1 == 4 && n2 == 3 && n3 == 2 && n4 == 2) return (14739259392 * sqrt(3) * Z) / 6131066257801.0;
    else if (n1 == 4 && n2 == 3 && n3 == 2 && n4 == 3) return (36645380390912 * sqrt(2) * Z) / 24984212408264457.0;
    else if (n1 == 4 && n2 == 3 && n3 == 2 && n4 == 4) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 4 && n2 == 3 && n3 == 3 && n4 == 1) return (30254432256 * Z) / 41426511213649.0;
    else if (n1 == 4 && n2 == 3 && n3 == 3 && n4 == 2) return (69158928384 * sqrt(2) * Z) / 34271896307633.0;
    else if (n1 == 4 && n2 == 3 && n3 == 3 && n4 == 3) return (2486755845603328 * Z) / (158298797548828125 * sqrt(3));
    else if (n1 == 4 && n2 == 3 && n3 == 3 && n4 == 4) return (2560158144 * Z) / 678223072849.0;
    else if (n1 == 4 && n2 == 3 && n3 == 4 && n4 == 1) return (1243165779 * sqrt(3) * Z) / 4564986729776.0;
    else if (n1 == 4 && n2 == 3 && n3 == 4 && n4 == 2) return (2496169683 * sqrt(1.5) * Z) / 1677721600000.0;
    else if (n1 == 4 && n2 == 3 && n3 == 4 && n4 == 3) return (621550729 * Z) / 13841287201.0;
    else if (n1 == 4 && n2 == 3 && n3 == 4 && n4 == 4) return (413631006610176000 * sqrt(3) * Z) / 249430673908303812379.0;
    else if (n1 == 4 && n2 == 4 && n3 == 1 && n4 == 1) return (22848 * Z) / 9765625.0;
    else if (n1 == 4 && n2 == 4 && n3 == 1 && n4 == 2) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 4 && n2 == 4 && n3 == 1 && n4 == 3) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 4 && n2 == 4 && n3 == 1 && n4 == 4) return (1804351488 * Z) / 6179146071875.0;
    else if (n1 == 4 && n2 == 4 && n3 == 2 && n4 == 1) return (39 * Z) / (32768 * sqrt(2));
    else if (n1 == 4 && n2 == 4 && n3 == 2 && n4 == 2) return (424 * Z) / 177147.0;
    else if (n1 == 4 && n2 == 4 && n3 == 2 && n4 == 3) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 4 && n2 == 4 && n3 == 2 && n4 == 4) return (21252608 * sqrt(2) * Z) / 35595703125.0;
    else if (n1 == 4 && n2 == 4 && n3 == 3 && n4 == 1) return (74450880 * sqrt(3) * Z) / 285311670611.0;
    else if (n1 == 4 && n2 == 4 && n3 == 3 && n4 == 2) return (145503 * sqrt(1.5) * Z) / 134217728.0;
    else if (n1 == 4 && n2 == 4 && n3 == 3 && n4 == 3) return (2560158144 * Z) / 678223072849.0;
    else if (n1 == 4 && n2 == 4 && n3 == 3 && n4 == 4) return (413631006610176000 * sqrt(3) * Z) / 249430673908303812379.0;
    else if (n1 == 4 && n2 == 4 && n3 == 4 && n4 == 1) return (1804351488 * Z) / 6179146071875.0;
    else if (n1 == 4 && n2 == 4 && n3 == 4 && n4 == 2) return (21252608 * sqrt(2) * Z) / 35595703125.0;
    else if (n1 == 4 && n2 == 4 && n3 == 4 && n4 == 3) return (413631006610176000 * sqrt(3) * Z) / 249430673908303812379.0;
    else if (n1 == 4 && n2 == 4 && n3 == 4 && n4 == 4) return (19541 * Z) / 524288;
    //    return N(0) * N(1) * N(2) * N(3) * I;











}

double hydrogenicOrbitals::get_sp_energy(int qnum) const {
    int n = qnums(qnum, 0);

    return -Z * Z / (2.0 * n * n);
}

void hydrogenicOrbitals::get_qnums() {
    int n2 = 5;
    qnums.set_size(n2, dim);

    int i = 0;
    int n = 1;
    while (i < n2) {

        for (int l = 0; l < n; l++) {
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



