#ifndef PRIMITIVEGTO_H
#define PRIMITIVEGTO_H

#include <iostream>
#include <armadillo>

class PrimitiveGTO
{
public:
    PrimitiveGTO(const double &exponent, const double &weight, const arma::rowvec &powers);

    double exponent() const;
    void setExponent(const double &exponent);

    double weight() const;
    void setWeight(const double &weight);

    arma::rowvec powers() const;
    void setPowers(const arma::rowvec &powers);

private:
    double m_exponent;
    double m_weight;
    arma::rowvec m_powers;

};
#endif // PRIMITIVEGTO_H
