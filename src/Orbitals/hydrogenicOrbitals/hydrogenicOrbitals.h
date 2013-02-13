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

    virtual void set_qnum_indie_terms(const Walker * walker, int i);

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

    virtual double get_variational_derivative(const Walker* walker, int n);

    void get_qnums();
    double L(int n, int l, double r) const;

    virtual double get_parameter(int n) {
        return *alpha;
    }

    virtual void set_parameter(double parameter, int n) {
        *alpha = parameter;
        *k = parameter*Z;
        *k2 = (*k)*(*k);
    }

};

#endif	/* HYDROGENICORBITALS_H */