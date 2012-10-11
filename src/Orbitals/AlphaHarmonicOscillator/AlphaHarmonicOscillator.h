/* 
 * File:   AlphaHarmonicOscillator.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 17:41
 */

#ifndef ALPHAHARMONICOSCILLATOR_H
#define	ALPHAHARMONICOSCILLATOR_H

class AlphaHarmonicOscillator : public Orbitals {
private:
    double *alpha;
    double *k;
    double *k2;
    double *w_over_a;
    
    arma::Mat<int> qnums;
    
    double w;

    virtual double get_parameter(int n);
    virtual void set_parameter(double parameter, int n);
    virtual double get_variational_derivative(const Walker* walker, int n) const;
    
    void get_qnums();
    double H(int n, double x) const;
    
public:
    AlphaHarmonicOscillator(GeneralParams &, VariationalParams &);
    AlphaHarmonicOscillator(GeneralParams &);

    friend class ExpandedBasis;
};

#endif	/* ALPHAHARMONICOSCILLATOR_H */