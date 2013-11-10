#include "h_321g.h"

H_321G::H_321G()
{

    rowvec exp1 = {5.4471780, 0.8245470};
    rowvec exp2 = {0.1831920};

    rowvec weight1 = {0.1562850, 0.9046910};
    rowvec weight2 = {1.0};


    weight1 = pow(2*exp1 / M_PI, 0.75)% weight1;
    weight2 = pow(2*exp2 / M_PI, 0.75)% weight2;

    rowvec powers = {0,0,0};

    PrimitiveGTO primitiveGTO_11(exp1(0),weight1(0) ,powers);
    PrimitiveGTO primitiveGTO_12(exp1(1),weight1(1) ,powers);
    PrimitiveGTO primitiveGTO_2(exp2(0), weight2(0)  ,powers);

    m_contractedGTO_1.addPrimitive(primitiveGTO_11);
    m_contractedGTO_1.addPrimitive(primitiveGTO_12);
    m_contractedGTO_2.addPrimitive(primitiveGTO_2);


    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);

    m_angularMomentum = 0;

}


int H_321G::getAngularMomentum() const
{
    return m_angularMomentum;
}
