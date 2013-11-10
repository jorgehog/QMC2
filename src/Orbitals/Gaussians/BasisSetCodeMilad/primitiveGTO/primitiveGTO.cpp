#include "primitiveGTO.h"

PrimitiveGTO::PrimitiveGTO(const double &exponent, const double &weight, const arma::rowvec &powers)
{
    setExponent(exponent);
    setWeight(weight);
    setPowers(powers);
}

double PrimitiveGTO::exponent() const
{
    return m_exponent;
}

void PrimitiveGTO::setExponent(const double &exponent)
{
    m_exponent = exponent;
}

double PrimitiveGTO::weight() const
{
    return m_weight;
}

void PrimitiveGTO::setWeight(const double &weight)
{
    m_weight = weight;
}

arma::rowvec PrimitiveGTO::powers() const
{
    return m_powers;
}

void PrimitiveGTO::setPowers(const arma::rowvec &powers)
{
    m_powers = powers;
}




