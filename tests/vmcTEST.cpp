/* 
 * File:   newsimpletest1.cpp
 * Author: jorgehog
 *
 * Created on 22.jun.2012, 12:52:55
 */
#include "../src/QMCheaders.h"

/*
 * Simple C++ Test Suite
 */

void test_nocol_nojast() {
    int n_p, dim, n_c;
    double alpha, w, dt, h;
    long random_seed;
    double var, E;

    dim = 2;
    n_c = 10000;
    alpha = 1;

    arma::vec NP(3);
    arma::vec W(2);
    W(0) = 1.0;
    W(1) = 0.5;
    NP(0) = 2;
    NP(1) = 6;
    NP(2) = 12;

    Orbitals* HO_basis;

    Potential* HO_pot;

    Jastrow* jastrow;

    Kinetics* kinetics;

    System* Quantum_Dots;

    Sampling* sample_method;

    VMC* vmc;

    dt = 0.5;

    h = 0.0001;

    random_seed = -time(NULL);


    std::cout << "Numerical BF no col no jast" << std::endl;
    bool success = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            n_p = NP(i);
            w = W(j);

            HO_basis = new oscillator_basis(n_p, dim, alpha, w);

            HO_pot = new Harmonic_osc(n_p, dim, w);

            jastrow = new No_Jastrow();

            kinetics = new Numerical(n_p, dim, h);

            Quantum_Dots = new Fermions(n_p, dim, HO_basis);

            Quantum_Dots->add_potential(HO_pot);

            sample_method = new Brute_Force(n_p, dim, dt, random_seed);

            vmc = new VMC(n_p, dim, n_c, jastrow, sample_method, Quantum_Dots, kinetics);


            vmc->run_method();

            var = vmc->get_var();
            E = vmc->get_energy();

            double E_CF = 1. / 3 * (i + 1)*(i + 2)*(2 * i + 3) * w;

            if (((E - E_CF)*(E - E_CF) < 0.1) && (var < 10E-5)) {
                std::cout << "n_p=" << n_p << " w=" << w << " success" << std::endl;
                std::cout << E << "\t" << var << std::endl;
                std::cout << endl;
            } else {
                std::cout << "n_p=" << n_p << " w=" << w << "failed" << std::endl;
                std::cout << E << "\t" << var << std::endl;
                std::cout << endl;
                success = false;
            }

        }

    }


    if (success == false) {
        std::cout << "%TEST_FAILED% time=0 testname=no_col_no_jast (nocolnojast_test) message=test failed" << std::endl;
    } else {
        std::cout << "%TEST_PASSED% time=0 testname=no_col_no_jsat (nocolnojast_test) message=yay" << std::endl;
    }
}

void test_col_jast_num() {
    int n_p, dim, n_c;
    double alpha, beta, w, dt, h;
    long random_seed;
    double var, E;

    dim = 2;
    n_c = 10000;

    arma::vec NP(3);
    arma::vec W(2);
    arma::mat E_ref = arma::zeros(3, 2);
    E_ref(0, 0) = 3.0;
    E_ref(0, 1) = 1.66;
    E_ref(1, 0) = 20.19;
    E_ref(1, 1) = 11.81;
    E_ref(2, 0) = 65.79;
    E_ref(2, 1) = 39.23;

    W(0) = 1.0;
    W(1) = 0.5;
    NP(0) = 2;
    NP(1) = 6;
    NP(2) = 12;

    Orbitals* HO_basis;

    Potential* HO_pot;

    Jastrow* jastrow;

    Kinetics* kinetics;

    System* Quantum_Dots;

    Sampling* sample_method;

    VMC* vmc;


    h = 0.0001;

    random_seed = -time(NULL);


    std::cout << "Numerical BF full" << std::endl;
    bool success = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {

            n_p = NP(i);
            w = W(j);

            initVMC(n_p, dim, w, dt, "QDots", "BF", alpha, beta);

            HO_basis = new oscillator_basis(n_p, dim, alpha, w);

            HO_pot = new Harmonic_osc(n_p, dim, w);
            jastrow = new Pade_Jastrow(n_p, dim, beta);

            kinetics = new Numerical(n_p, dim, h);

            Quantum_Dots = new Fermions(n_p, dim, HO_basis);
            Quantum_Dots->add_potential(HO_pot);
            Quantum_Dots->add_potential(new Coulomb(n_p, dim));

            sample_method = new Brute_Force(n_p, dim, dt, random_seed);

            vmc = new VMC(n_p, dim, n_c, jastrow, sample_method, Quantum_Dots, kinetics);

            vmc->run_method();

            var = vmc->get_var();
            E = vmc->get_energy();

            if ((E - E_ref(i, j))*(E - E_ref(i, j)) < 0.5 * n_p) {
                std::cout << "n_p=" << n_p << " w=" << w << " success" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
            } else {
                std::cout << "n_p=" << n_p << " w=" << w << "failed" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
                success = false;
            }

        }

    }

    if (success == false) {
        std::cout << "%TEST_FAILED% time=0 testname=test_col_jast_num (coljastnum_test) message=test failed" << std::endl;
    }
}

void test_col_jast_CF() {
    int n_p, dim, n_c;
    double alpha, beta, w, dt, h;
    long random_seed;
    double var, E;

    dim = 2;
    n_c = 10000;

    arma::vec NP(3);
    arma::vec W(2);
    arma::mat E_ref = arma::zeros(3, 2);
    E_ref(0, 0) = 3.0;
    E_ref(0, 1) = 1.66;
    E_ref(1, 0) = 20.19;
    E_ref(1, 1) = 11.81;
    E_ref(2, 0) = 65.79;
    E_ref(2, 1) = 39.23;

    W(0) = 1.0;
    W(1) = 0.5;
    NP(0) = 2;
    NP(1) = 6;
    NP(2) = 12;

    Orbitals* HO_basis;

    Potential* HO_pot;

    Jastrow* jastrow;

    Kinetics* kinetics;

    System* Quantum_Dots;

    Sampling* sample_method;

    VMC* vmc;


    random_seed = -time(NULL);


    std::cout << "Closed Form BF full" << std::endl;
    bool success = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {

            n_p = NP(i);
            w = W(j);

            initVMC(n_p, dim, w, dt, "QDots", "BF", alpha, beta);

            HO_basis = new oscillator_basis(n_p, dim, alpha, w);

            HO_pot = new Harmonic_osc(n_p, dim, w);

            jastrow = new Pade_Jastrow(n_p, dim, beta);

            kinetics = new Closed_form(n_p, dim);

            Quantum_Dots = new Fermions(n_p, dim, HO_basis);
            Quantum_Dots->add_potential(HO_pot);
            Quantum_Dots->add_potential(new Coulomb(n_p, dim));

            sample_method = new Brute_Force(n_p, dim, dt, random_seed);

            vmc = new VMC(n_p, dim, n_c, jastrow, sample_method, Quantum_Dots, kinetics);

            vmc->run_method();

            var = vmc->get_var();
            E = vmc->get_energy();

            if ((E - E_ref(i, j))*(E - E_ref(i, j)) < 0.5 * n_p) {
                std::cout << "n_p=" << n_p << " w=" << w << " success" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
            } else {
                std::cout << "n_p=" << n_p << " w=" << w << "failed" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
                success = false;
            }

        }
    }

    if (success == false) {
        std::cout << "%TEST_FAILED% time=0 testname=test_col_jast_CF (coljastCF_test) message=test failed" << std::endl;
    }
}

void test_ISnum_nocoljast() {
    int n_p, dim, n_c;
    double alpha, w, dt, h;
    long random_seed;
    double var, E;

    dim = 2;
    n_c = 10000;
    alpha = 1;

    arma::vec NP(3);
    arma::vec W(2);
    W(0) = 1.0;
    W(1) = 0.5;
    NP(0) = 2;
    NP(1) = 6;
    NP(2) = 12;

    Orbitals* HO_basis;

    Potential* HO_pot;

    Jastrow* jastrow;

    Kinetics* kinetics;

    System* Quantum_Dots;

    Sampling* sample_method;

    VMC* vmc;

    dt = 0.5;

    h = 0.0001;

    random_seed = -time(NULL);


    std::cout << "Numerical IS no col no jast" << std::endl;
    bool success = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            n_p = NP(i);
            w = W(j);

            HO_basis = new oscillator_basis(n_p, dim, alpha, w);

            HO_pot = new Harmonic_osc(n_p, dim, w);

            jastrow = new No_Jastrow();

            kinetics = new Numerical(n_p, dim, h);

            Quantum_Dots = new Fermions(n_p, dim, HO_basis);
            Quantum_Dots->add_potential(HO_pot);

            sample_method = new Importance(n_p, dim, dt, random_seed);

            vmc = new VMC(n_p, dim, n_c, jastrow, sample_method, Quantum_Dots, kinetics);
            
            cout << "foer" << endl;
            vmc->run_method();
            cout << "etter" << endl;
            var = vmc->get_var();
            E = vmc->get_energy();
            double E_CF = 1. / 3 * (i + 1)*(i + 2)*(2 * i + 3) * w;

            if (((E - E_CF)*(E - E_CF) < 0.1) && (var < 10E-5)) {
                std::cout << "n_p=" << n_p << " w=" << w << " success" << std::endl;
                std::cout << E << "\t" << var << std::endl;
                std::cout << endl;
            } else {
                std::cout << "n_p=" << n_p << " w=" << w << "failed" << std::endl;
                std::cout << E << "\t" << var << std::endl;
                std::cout << endl;
                success = false;
            }

        }

    }


    if (success == false) {
        std::cout << "%TEST_FAILED% time=0 testname=test_ISnum_nocoljast (ISnum_nocoljast_test) message=test failed" << std::endl;
    }
}

void test_ISnum() {
    int n_p, dim, n_c;
    double alpha, beta, w, dt, h;
    long random_seed;
    double var, E;

    dim = 2;
    n_c = 10000;

    arma::vec NP(3);
    arma::vec W(2);
    arma::mat E_ref = arma::zeros(3, 2);
    E_ref(0, 0) = 3.0;
    E_ref(0, 1) = 1.66;
    E_ref(1, 0) = 20.19;
    E_ref(1, 1) = 11.81;
    E_ref(2, 0) = 65.79;
    E_ref(2, 1) = 39.23;

    W(0) = 1.0;
    W(1) = 0.5;
    NP(0) = 2;
    NP(1) = 6;
    NP(2) = 12;

    Orbitals* HO_basis;

    Potential* HO_pot;

    Jastrow* jastrow;

    Kinetics* kinetics;

    System* Quantum_Dots;

    Sampling* sample_method;

    VMC* vmc;


    random_seed = -time(NULL);


    std::cout << "Numerical IS full" << std::endl;
    bool success = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {

            n_p = NP(i);
            w = W(j);

            initVMC(n_p, dim, w, dt, "QDots", "IS", alpha, beta);

            HO_basis = new oscillator_basis(n_p, dim, alpha, w);

            HO_pot = new Harmonic_osc(n_p, dim, w);

            jastrow = new Pade_Jastrow(n_p, dim, beta);

            kinetics = new Numerical(n_p, dim, 0.001);

            Quantum_Dots = new Fermions(n_p, dim, HO_basis);
            Quantum_Dots->add_potential(HO_pot);
            Quantum_Dots->add_potential(new Coulomb(n_p, dim));

            sample_method = new Importance(n_p, dim, dt, random_seed);

            vmc = new VMC(n_p, dim, n_c, jastrow, sample_method, Quantum_Dots, kinetics);

            vmc->run_method();

            var = vmc->get_var();
            E = vmc->get_energy();

            if ((E - E_ref(i, j))*(E - E_ref(i, j)) < 0.5 * n_p) {
                std::cout << "n_p=" << n_p << " w=" << w << " success" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
            } else {
                std::cout << "n_p=" << n_p << " w=" << w << "failed" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
                success = false;
            }

        }
    }

    if (success == false) {
        std::cout << "%TEST_FAILED% time=0 testname=test_ISnum (ISnum_test) message=test failed" << std::endl;
    }
}

void test_ISCF() {
    int n_p, dim, n_c;
    double alpha, beta, w, dt, h;
    long random_seed;
    double var, E;

    dim = 2;
    n_c = 10000;

    arma::vec NP(3);
    arma::vec W(2);
    arma::mat E_ref = arma::zeros(3, 2);
    E_ref(0, 0) = 3.0;
    E_ref(0, 1) = 1.66;
    E_ref(1, 0) = 20.19;
    E_ref(1, 1) = 11.81;
    E_ref(2, 0) = 65.79;
    E_ref(2, 1) = 39.23;

    W(0) = 1.0;
    W(1) = 0.5;
    NP(0) = 2;
    NP(1) = 6;
    NP(2) = 12;

    Orbitals* HO_basis;

    Potential* HO_pot;

    Jastrow* jastrow;

    Kinetics* kinetics;

    System* Quantum_Dots;

    Sampling* sample_method;

    VMC* vmc;


    random_seed = -time(NULL);


    std::cout << "CF IS full" << std::endl;
    bool success = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {

            n_p = NP(i);
            w = W(j);

            initVMC(n_p, dim, w, dt, "QDots", "IS", alpha, beta);

            HO_basis = new oscillator_basis(n_p, dim, alpha, w);

            HO_pot = new Harmonic_osc(n_p, dim, w);

            jastrow = new Pade_Jastrow(n_p, dim, beta);

            kinetics = new Closed_form(n_p, dim);

            Quantum_Dots = new Fermions(n_p, dim, HO_basis);
            Quantum_Dots->add_potential(HO_pot);
            Quantum_Dots->add_potential(new Coulomb(n_p, dim));

            sample_method = new Importance(n_p, dim, dt, random_seed);

            vmc = new VMC(n_p, dim, n_c, jastrow, sample_method, Quantum_Dots, kinetics);

            vmc->run_method();

            var = vmc->get_var();
            E = vmc->get_energy();

            if ((E - E_ref(i, j))*(E - E_ref(i, j)) < 0.5 * n_p) {
                std::cout << "n_p=" << n_p << " w=" << w << " success" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
            } else {
                std::cout << "n_p=" << n_p << " w=" << w << "failed" << std::endl;
                std::cout << E << "\t" << E_ref(i, j) << std::endl;
                std::cout << endl;
                success = false;
            }

        }
    }

    if (success == false) {
        std::cout << "%TEST_FAILED% time=0 testname=test_ISCF (ISCF_test) message=test failed" << std::endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% vmcTEST" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    /*
     
     
     
     */

    std::cout << "%TEST_STARTED% no_col_no_jast (nocolnojast_test)" << std::endl;
    test_nocol_nojast();
    std::cout << "%TEST_FINISHED% time=0 no_col_no_jast (nocolnojast_test)" << std::endl;

    std::cout << "%TEST_STARTED% test_col_jast_num (coljastnum_test)" << std::endl;
    test_col_jast_num();
    std::cout << "%TEST_FINISHED% time=0 test_col_jast_num (coljastnum_test)" << std::endl;

    std::cout << "%TEST_STARTED% test_col_jast_CF (coljastCF_test)" << std::endl;
    test_col_jast_CF();
    std::cout << "%TEST_FINISHED% time=0 test_col_jast_CF (coljastCF_test)" << std::endl;

    std::cout << "%TEST_STARTED% test_ISnum_nocoljast (ISnum_nocoljast_test)" << std::endl;
    test_ISnum_nocoljast();
    std::cout << "%TEST_FINISHED% time=0 test_ISnumnocoljast (ISnum_nocoljast_test)" << std::endl;

    std::cout << "%TEST_STARTED% test_ISnum (ISnum_test)" << std::endl;
    test_ISnum();
    std::cout << "%TEST_FINISHED% time=0 test_ISnum (ISnum_test)" << std::endl;

    std::cout << "%TEST_STARTED% test_ISCF (ISCF_test)" << std::endl;
    test_ISCF();
    std::cout << "%TEST_FINISHED% time=0 test_ISCF (ISCF_test)" << std::endl;
    /*
     
     
     */


    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

