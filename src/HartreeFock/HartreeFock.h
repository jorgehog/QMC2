/* 
 * File:   HartreeFock.h
 * Author: jorgmeister
 *
 * Created on April 4, 2013, 2:46 PM
 */

#ifndef HARTREEFOCK_H
#define	HARTREEFOCK_H

#include <armadillo>
class Orbitals;

class HartreeFock {
public:
    HartreeFock(int m, Orbitals* sp_basis, double thresh = 1e-5);

    void run_method();

private:

    double thresh;

    int n_p;
    int m;
    int m2;
    
    Orbitals* sp_basis;

    arma::mat H_hf;

    arma::mat C;
    arma::vec e_HF;

    arma::vec e_sp;
    arma::field<arma::mat> V;


    void setup_V();
    void setup_sp_e();
    
    double get_V(int a, int b, int g, int d);

    void make_Hamiltonian();

    double calc_HF_E();
    
    void sort();

};

#endif	/* HARTREEFOCK_H */

