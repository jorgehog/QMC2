/* 
 * File:   AlphaHarmonicOscillator.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 17:41
 */

#ifndef ALPHAHARMONICOSCILLATOR_H
#define	ALPHAHARMONICOSCILLATOR_H

class AlphaHarmonicOscillator : public Orbitals {
public:
    AlphaHarmonicOscillator(GeneralParams &, VariationalParams &);
    AlphaHarmonicOscillator(GeneralParams &);

    virtual void set_qnum_indie_terms(const Walker * walker, int i) {
        *exp_factor = exp(-0.5 * (*k2) * walker->get_r_i2(i));
    }

    friend class ExpandedBasis;

protected:
    double *alpha;
    double *k;
    double *k2;
    double *w_over_a;

    double *exp_factor;

    arma::Mat<int> qnums;

    double w;

    virtual double get_variational_derivative(const Walker* walker, int n);

    void get_qnums();
    double H(int n, double x) const;

    virtual double get_parameter(int n) {
        return *alpha;
    }

    virtual void set_parameter(double parameter, int n) {
        *alpha = parameter;
        *k2 = parameter*w;
        *k = sqrt(*k2);
    }

};

#endif	/* ALPHAHARMONICOSCILLATOR_H */