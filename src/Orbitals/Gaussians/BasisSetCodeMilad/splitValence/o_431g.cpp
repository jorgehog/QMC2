#include "o_431g.h"

O_431G::O_431G()
{

    rowvec exp1 = {883.2728600, 133.1292800, 29.9064080, 7.9786772 };
    rowvec exp2 = {16.1944470, 3.7800860, 1.0709836};
    rowvec exp3 = {0.2838798 };
    rowvec exp4 = {16.1944470, 3.7800860, 1.0709836};
    rowvec exp5 = {0.2838798 };

    rowvec weight1 = {0.0175506, 0.1228292, 0.4348836, 0.5600108};
    rowvec weight2 = {-0.1134010, -0.1772865, 1.1504079 };
    rowvec weight3 = {1.0};
    rowvec weight4 = {0.0685453 ,0.3312254, 0.7346079};
    rowvec weight5 = {1.0};

    rowvec powers1 = {0, 0, 0};
    rowvec powers2 = {1, 0, 0};
    rowvec powers3 = {0, 1, 0};
    rowvec powers4 = {0, 0, 1};


    weight1 = pow(2*exp1 / M_PI, 0.75)% weight1;
    weight2 = pow(2*exp2 / M_PI, 0.75)% weight2;
    weight3 = pow(2*exp3 / M_PI, 0.75)% weight3;

    weight4 = pow(2*exp4 / M_PI, 0.75) * 2 % sqrt(exp4) % weight4;
    weight5 = pow(2*exp5 / M_PI, 0.75) * 2 % sqrt(exp5) % weight5;


    //1s
    PrimitiveGTO primitiveGTO_11(exp1(0),weight1(0) ,powers1);
    PrimitiveGTO primitiveGTO_12(exp1(1),weight1(1) ,powers1);
    PrimitiveGTO primitiveGTO_13(exp1(2),weight1(2) ,powers1);
    PrimitiveGTO primitiveGTO_14(exp1(3),weight1(3) ,powers1);
    m_contractedGTO_1.addPrimitive(primitiveGTO_11);
    m_contractedGTO_1.addPrimitive(primitiveGTO_12);
    m_contractedGTO_1.addPrimitive(primitiveGTO_13);
    m_contractedGTO_1.addPrimitive(primitiveGTO_14);


    //2s
    PrimitiveGTO primitiveGTO_21(exp2(0),weight2(0) ,powers1);
    PrimitiveGTO primitiveGTO_22(exp2(1),weight2(1) ,powers1);
    PrimitiveGTO primitiveGTO_23(exp2(2),weight2(2) ,powers1);

    PrimitiveGTO primitiveGTO_3(exp3(0),weight3(0) ,powers1);

    m_contractedGTO_2.addPrimitive(primitiveGTO_21);
    m_contractedGTO_2.addPrimitive(primitiveGTO_22);
    m_contractedGTO_2.addPrimitive(primitiveGTO_23);

    m_contractedGTO_3.addPrimitive(primitiveGTO_3);

    //2px
    PrimitiveGTO primitiveGTO_41(exp4(0),weight4(0) ,powers2);
    PrimitiveGTO primitiveGTO_42(exp4(1),weight4(1) ,powers2);
    PrimitiveGTO primitiveGTO_43(exp4(2),weight4(2) ,powers2);

    PrimitiveGTO primitiveGTO_5(exp5(0),weight5(0) ,powers2);

    m_contractedGTO_4.addPrimitive(primitiveGTO_41);
    m_contractedGTO_4.addPrimitive(primitiveGTO_42);
    m_contractedGTO_4.addPrimitive(primitiveGTO_43);

    m_contractedGTO_5.addPrimitive(primitiveGTO_5);

    //2py
    PrimitiveGTO primitiveGTO_61(exp4(0),weight4(0) ,powers3);
    PrimitiveGTO primitiveGTO_62(exp4(1),weight4(1) ,powers3);
    PrimitiveGTO primitiveGTO_63(exp4(2),weight4(2) ,powers3);

    PrimitiveGTO primitiveGTO_7(exp5(0),weight5(0) ,powers3);

    m_contractedGTO_6.addPrimitive(primitiveGTO_61);
    m_contractedGTO_6.addPrimitive(primitiveGTO_62);
    m_contractedGTO_6.addPrimitive(primitiveGTO_63);

    m_contractedGTO_7.addPrimitive(primitiveGTO_7);

    //2pz
    PrimitiveGTO primitiveGTO_81(exp4(0),weight4(0) ,powers4);
    PrimitiveGTO primitiveGTO_82(exp4(1),weight4(1) ,powers4);
    PrimitiveGTO primitiveGTO_83(exp4(2),weight4(2) ,powers4);

    PrimitiveGTO primitiveGTO_9(exp5(0),weight5(0) ,powers4);

    m_contractedGTO_8.addPrimitive(primitiveGTO_81);
    m_contractedGTO_8.addPrimitive(primitiveGTO_82);
    m_contractedGTO_8.addPrimitive(primitiveGTO_83);

    m_contractedGTO_9.addPrimitive(primitiveGTO_9);


    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);
    m_contractedGTOs.push_back(m_contractedGTO_3);
    m_contractedGTOs.push_back(m_contractedGTO_4);
    m_contractedGTOs.push_back(m_contractedGTO_5);
    m_contractedGTOs.push_back(m_contractedGTO_6);
    m_contractedGTOs.push_back(m_contractedGTO_7);
    m_contractedGTOs.push_back(m_contractedGTO_8);
    m_contractedGTOs.push_back(m_contractedGTO_9);


    m_angularMomentum = 1;

}


int O_431G::getAngularMomentum() const
{
    return m_angularMomentum;
}
