#include "h_431g.h"

H_431G::H_431G()
{


    rowvec exp1 = {18.7311370, 2.8253944, 0.6401217};
    rowvec exp2 = {0.1612778 };

    rowvec weight1 = {0.0334946, 0.2347269, 0.8137573};
    rowvec weight2 = {1.0};


    weight1 = pow(2*exp1 / M_PI, 0.75)% weight1;
    weight2 = pow(2*exp2 / M_PI, 0.75)% weight2;

    rowvec powers = {0,0,0};

    PrimitiveGTO primitiveGTO_11(exp1(0),weight1(0) ,powers);
    PrimitiveGTO primitiveGTO_12(exp1(1),weight1(1) ,powers);
    PrimitiveGTO primitiveGTO_13(exp1(2),weight1(2) ,powers);

    PrimitiveGTO primitiveGTO_2(exp2(0), weight2(0) ,powers);

    m_contractedGTO_1.addPrimitive(primitiveGTO_11);
    m_contractedGTO_1.addPrimitive(primitiveGTO_12);
    m_contractedGTO_1.addPrimitive(primitiveGTO_13);
    m_contractedGTO_2.addPrimitive(primitiveGTO_2);


    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);

    m_angularMomentum = 0;
}


int H_431G::getAngularMomentum() const
{
    return m_angularMomentum;
}
