/* 
 * File:   alphaHO.h
 * Author: jorgehog
 *
 * Created on 27. juni 2012, 13:17
 */

#ifndef ALPHAHO_H
#define ALPHAHO_H

class alphaHO : public BasisFunctions {
protected:

    double* k;
    double* k2;
    double* exp_factor;

    double H;
    double x2;
    double y2;

public:

    alphaHO(double* k, double* k2, double* exp_factor);

};

//START OF GENERATED FUNCTIONS

class alphaHO_0 : public alphaHO {
public:

    alphaHO_0(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_0_x : public alphaHO {
public:

    dell_alphaHO_0_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_0_y : public alphaHO {
public:

    dell_alphaHO_0_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_0 : public alphaHO {
public:

    lapl_alphaHO_0(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_1 : public alphaHO {
public:

    alphaHO_1(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_1_x : public alphaHO {
public:

    dell_alphaHO_1_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_1_y : public alphaHO {
public:

    dell_alphaHO_1_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_1 : public alphaHO {
public:

    lapl_alphaHO_1(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_2 : public alphaHO {
public:

    alphaHO_2(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_2_x : public alphaHO {
public:

    dell_alphaHO_2_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_2_y : public alphaHO {
public:

    dell_alphaHO_2_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_2 : public alphaHO {
public:

    lapl_alphaHO_2(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_3 : public alphaHO {
public:

    alphaHO_3(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_3_x : public alphaHO {
public:

    dell_alphaHO_3_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_3_y : public alphaHO {
public:

    dell_alphaHO_3_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_3 : public alphaHO {
public:

    lapl_alphaHO_3(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_4 : public alphaHO {
public:

    alphaHO_4(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_4_x : public alphaHO {
public:

    dell_alphaHO_4_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_4_y : public alphaHO {
public:

    dell_alphaHO_4_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_4 : public alphaHO {
public:

    lapl_alphaHO_4(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_5 : public alphaHO {
public:

    alphaHO_5(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_5_x : public alphaHO {
public:

    dell_alphaHO_5_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_5_y : public alphaHO {
public:

    dell_alphaHO_5_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_5 : public alphaHO {
public:

    lapl_alphaHO_5(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_6 : public alphaHO {
public:

    alphaHO_6(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_6_x : public alphaHO {
public:

    dell_alphaHO_6_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_6_y : public alphaHO {
public:

    dell_alphaHO_6_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_6 : public alphaHO {
public:

    lapl_alphaHO_6(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_7 : public alphaHO {
public:

    alphaHO_7(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_7_x : public alphaHO {
public:

    dell_alphaHO_7_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_7_y : public alphaHO {
public:

    dell_alphaHO_7_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_7 : public alphaHO {
public:

    lapl_alphaHO_7(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_8 : public alphaHO {
public:

    alphaHO_8(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_8_x : public alphaHO {
public:

    dell_alphaHO_8_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_8_y : public alphaHO {
public:

    dell_alphaHO_8_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_8 : public alphaHO {
public:

    lapl_alphaHO_8(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_9 : public alphaHO {
public:

    alphaHO_9(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_9_x : public alphaHO {
public:

    dell_alphaHO_9_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_9_y : public alphaHO {
public:

    dell_alphaHO_9_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_9 : public alphaHO {
public:

    lapl_alphaHO_9(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_10 : public alphaHO {
public:

    alphaHO_10(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_10_x : public alphaHO {
public:

    dell_alphaHO_10_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_10_y : public alphaHO {
public:

    dell_alphaHO_10_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_10 : public alphaHO {
public:

    lapl_alphaHO_10(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_11 : public alphaHO {
public:

    alphaHO_11(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_11_x : public alphaHO {
public:

    dell_alphaHO_11_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_11_y : public alphaHO {
public:

    dell_alphaHO_11_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_11 : public alphaHO {
public:

    lapl_alphaHO_11(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_12 : public alphaHO {
public:

    alphaHO_12(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_12_x : public alphaHO {
public:

    dell_alphaHO_12_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_12_y : public alphaHO {
public:

    dell_alphaHO_12_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_12 : public alphaHO {
public:

    lapl_alphaHO_12(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_13 : public alphaHO {
public:

    alphaHO_13(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_13_x : public alphaHO {
public:

    dell_alphaHO_13_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_13_y : public alphaHO {
public:

    dell_alphaHO_13_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_13 : public alphaHO {
public:

    lapl_alphaHO_13(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

class alphaHO_14 : public alphaHO {
public:

    alphaHO_14(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_14_x : public alphaHO {
public:

    dell_alphaHO_14_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_alphaHO_14_y : public alphaHO {
public:

    dell_alphaHO_14_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_alphaHO_14 : public alphaHO {
public:

    lapl_alphaHO_14(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

//-----------------------------------//

//END OF GENERATED FUNCTIONS

#endif	/* ALPHAHO_H */
