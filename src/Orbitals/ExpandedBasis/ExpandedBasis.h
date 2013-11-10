/* 
 * File:   ExpandedBasis.h
 * Author: jorgmeister
 *
 * Created on October 9, 2012, 5:09 PM
 */

#ifndef EXPANDEDBASIS_H
#define	EXPANDEDBASIS_H

#include "../Orbitals.h"

#include "../../structs.h"
#include "../../Walker/Walker.h"

class ExpandedBasis : public Orbitals {
public:
    ExpandedBasis(GeneralParams & gp) :
        Orbitals(gp.n_p, gp.dim) {

    }

    void setBasis(Orbitals* basis) {
        name = "Expanded" + basis->getName();
        this->basis = basis;
    }

    void setCoeffs(const arma::mat &coeffs) {
        this->coeffs = coeffs;
        this->basis_size = coeffs.n_cols;
    }

    arma::mat getCoeffs() {
        return coeffs;
    }

    void set_qnum_indie_terms(Walker* walker, int i) {
        basis->set_qnum_indie_terms(walker, i);
    }


    double phi(const Walker* walker, int particle, int q_num) {
        double value = 0;

        //Dividing basis_size by half assuming a two-level system.
        //In case of Bosons, expanding s.p. w.f. does not make sence.
        for (int m = 0; m < basis_size; m++) {
            value += coeffs(q_num, m) * basis->phi(walker, particle, m);
        }
        return value;
    }

    double del_phi(const Walker* walker, int particle, int q_num, int d) {
        double value = 0;
        for (int m = 0; m < basis_size; m++) {
            value += coeffs(particle, m) * basis->del_phi(walker, particle, q_num, d);
        }

        return value;
    }

    double lapl_phi(const Walker* walker, int particle, int q_num) {
        double value = 0;
        for (int m = 0; m < basis_size; m++) {
            value += coeffs(particle, m) * basis->lapl_phi(walker, particle, q_num);
        }

        return value;
    }

protected:
    int basis_size;

    arma::mat coeffs;

    Orbitals* basis = NULL;

};

#endif	/* EXPANDEDBASIS_H */

