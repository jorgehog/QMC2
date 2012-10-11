/* 
 * File:   AlphaHarmonicOscillatorOld.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 17:41
 */

#ifndef ALPHAHARMONICOSCILLATOROLD_H
#define	ALPHAHARMONICOSCILLATOROLD_H

class AlphaHarmonicOscillatorOld : public Orbitals {
private:
    double *alpha;
    
    double w;

    virtual double get_parameter(int n);
    virtual void set_parameter(double parameter, int n);
    virtual double get_variational_derivative(const Walker* walker, int n) const;
    
public:
    AlphaHarmonicOscillatorOld(GeneralParams &, VariationalParams &);

    friend class ExpandedBasis;
};

#endif	/* ALPHAHARMONICOSCILLATOROLD_H */