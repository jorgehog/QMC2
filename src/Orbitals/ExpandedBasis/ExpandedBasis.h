#pragma once


#include "../Orbitals.h"

#include "../../structs.h"
#include "../../Walker/Walker.h"


namespace QMC2
{


class ExpandedBasis : public Orbitals {
public:

    ExpandedBasis() {}

    ExpandedBasis(Orbitals* basis, const arma::mat &coeffs)
    {
        setBasis(basis);
        setCoeffs(coeffs);
    }

    void setBasis(Orbitals* basis) {
        name = "Expanded" + basis->getName();
        this->basis = basis;
    }

    void setCoeffs(const arma::mat &coeffs) {
        this->coeffs = coeffs;
        this->basis_size = coeffs.n_rows;
    }

    void set_parameter(double parameter, int n) {
        basis->set_parameter(parameter, n);
    }

    double get_parameter(int n){
        return basis->get_parameter(n);
    }

    const arma::mat & getCoeffs() {
        return coeffs;
    }

    void set_qnum_indie_terms(Walker* walker, int i) {
        basis->set_qnum_indie_terms(walker, i);
    }

    double get_dell_alpha_phi(Walker *walker, int p, int q_num, int n) {
        double value = 0;

        for (int m = 0; m < basis_size; m++) {
            value += coeffs(m, q_num) * basis->get_dell_alpha_phi(walker, p, q_num, n);
        }
        return value;
    }

    double phi(const Walker* walker, int particle, int q_num) {
        double value = 0;


        for (int m = 0; m < basis_size; m++) {
            value += coeffs(m, q_num) * basis->phi(walker, particle, m);
        }
        return value;
    }

    double del_phi(const Walker* walker, int particle, int q_num, int d) {
        double value = 0;
        for (int m = 0; m < basis_size; m++) {
            value += coeffs(m, q_num) * basis->del_phi(walker, particle, m, d);
        }

        return value;
    }

    double lapl_phi(const Walker* walker, int particle, int q_num) {
        double value = 0;
        for (int m = 0; m < basis_size; m++) {
            value += coeffs(m, q_num) * basis->lapl_phi(walker, particle, m);
        }

        return value;
    }

protected:
    int basis_size;

    arma::mat coeffs;

    Orbitals* basis;

};

}
