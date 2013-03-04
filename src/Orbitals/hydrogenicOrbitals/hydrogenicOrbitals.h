/* 
 * File:   hydrogenicOrbitals.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 17:41
 */

#ifndef HYDROGENICORBITALS_H
#define	HYDROGENICORBITALS_H

class hydrogenicOrbitals : public Orbitals {
public:
    hydrogenicOrbitals(GeneralParams &, VariationalParams &);
    hydrogenicOrbitals(GeneralParams &);

    void set_qnum_indie_terms(const Walker * walker, int i);

    friend class ExpandedBasis;

protected:
    double *alpha;
    double *k;
    double *k2;
    double *r22d;
    double *r2d;
    
    
    double *exp_factor_n1;
    double *exp_factor_n2;

    arma::Mat<int> qnums;

    int Z;

    double get_variational_derivative(const Walker* walker, int n);

    double get_dell_alpha_phi(const Walker* walker, int qnum, int i);

    double get_parameter(int n) {
        return *alpha;
    }

    void set_parameter(double parameter, int n) {
        *alpha = parameter;
        *k = parameter*Z;
        *k2 = (*k)*(*k);
    }

};

#endif	/* HYDROGENICORBITALS_H */