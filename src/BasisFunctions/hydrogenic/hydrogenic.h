#pragma once

#include "../BasisFunctions.h"

//Superclass 

class hydrogenic : public BasisFunctions {
protected:

    double* k;
    double* k2;
    double* exp_factor;
    double* r22d;
    double* r2d;

    double psi;
    double x;
    double y;
    double z;
    double x2;
    double y2;
    double z2;

public:

    hydrogenic(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);

};


/*
    Subclasses
*/


class hydrogenic_0 : public hydrogenic {
public:

    hydrogenic_0(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_0_x : public hydrogenic {
public:

    dell_hydrogenic_0_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_0_y : public hydrogenic {
public:

    dell_hydrogenic_0_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_0_z : public hydrogenic {
public:

    dell_hydrogenic_0_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_0 : public hydrogenic {
public:

    lapl_hydrogenic_0(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 0  -------------------------
*/

class hydrogenic_1 : public hydrogenic {
public:

    hydrogenic_1(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_1_x : public hydrogenic {
public:

    dell_hydrogenic_1_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_1_y : public hydrogenic {
public:

    dell_hydrogenic_1_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_1_z : public hydrogenic {
public:

    dell_hydrogenic_1_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_1 : public hydrogenic {
public:

    lapl_hydrogenic_1(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 1  -------------------------
*/

class hydrogenic_2 : public hydrogenic {
public:

    hydrogenic_2(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_2_x : public hydrogenic {
public:

    dell_hydrogenic_2_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_2_y : public hydrogenic {
public:

    dell_hydrogenic_2_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_2_z : public hydrogenic {
public:

    dell_hydrogenic_2_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_2 : public hydrogenic {
public:

    lapl_hydrogenic_2(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 2  -------------------------
*/

class hydrogenic_3 : public hydrogenic {
public:

    hydrogenic_3(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_3_x : public hydrogenic {
public:

    dell_hydrogenic_3_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_3_y : public hydrogenic {
public:

    dell_hydrogenic_3_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_3_z : public hydrogenic {
public:

    dell_hydrogenic_3_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_3 : public hydrogenic {
public:

    lapl_hydrogenic_3(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 3  -------------------------
*/

class hydrogenic_4 : public hydrogenic {
public:

    hydrogenic_4(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_4_x : public hydrogenic {
public:

    dell_hydrogenic_4_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_4_y : public hydrogenic {
public:

    dell_hydrogenic_4_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_4_z : public hydrogenic {
public:

    dell_hydrogenic_4_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_4 : public hydrogenic {
public:

    lapl_hydrogenic_4(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 4  -------------------------
*/

class hydrogenic_5 : public hydrogenic {
public:

    hydrogenic_5(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_5_x : public hydrogenic {
public:

    dell_hydrogenic_5_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_5_y : public hydrogenic {
public:

    dell_hydrogenic_5_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_5_z : public hydrogenic {
public:

    dell_hydrogenic_5_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_5 : public hydrogenic {
public:

    lapl_hydrogenic_5(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 5  -------------------------
*/

class hydrogenic_6 : public hydrogenic {
public:

    hydrogenic_6(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_6_x : public hydrogenic {
public:

    dell_hydrogenic_6_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_6_y : public hydrogenic {
public:

    dell_hydrogenic_6_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_6_z : public hydrogenic {
public:

    dell_hydrogenic_6_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_6 : public hydrogenic {
public:

    lapl_hydrogenic_6(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 6  -------------------------
*/

class hydrogenic_7 : public hydrogenic {
public:

    hydrogenic_7(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_7_x : public hydrogenic {
public:

    dell_hydrogenic_7_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_7_y : public hydrogenic {
public:

    dell_hydrogenic_7_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_7_z : public hydrogenic {
public:

    dell_hydrogenic_7_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_7 : public hydrogenic {
public:

    lapl_hydrogenic_7(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 7  -------------------------
*/

class hydrogenic_8 : public hydrogenic {
public:

    hydrogenic_8(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_8_x : public hydrogenic {
public:

    dell_hydrogenic_8_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_8_y : public hydrogenic {
public:

    dell_hydrogenic_8_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_8_z : public hydrogenic {
public:

    dell_hydrogenic_8_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_8 : public hydrogenic {
public:

    lapl_hydrogenic_8(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 8  -------------------------
*/

class hydrogenic_9 : public hydrogenic {
public:

    hydrogenic_9(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_9_x : public hydrogenic {
public:

    dell_hydrogenic_9_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_9_y : public hydrogenic {
public:

    dell_hydrogenic_9_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_9_z : public hydrogenic {
public:

    dell_hydrogenic_9_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_9 : public hydrogenic {
public:

    lapl_hydrogenic_9(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 9  -------------------------
*/

class hydrogenic_10 : public hydrogenic {
public:

    hydrogenic_10(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_10_x : public hydrogenic {
public:

    dell_hydrogenic_10_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_10_y : public hydrogenic {
public:

    dell_hydrogenic_10_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_10_z : public hydrogenic {
public:

    dell_hydrogenic_10_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_10 : public hydrogenic {
public:

    lapl_hydrogenic_10(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 10  -------------------------
*/

class hydrogenic_11 : public hydrogenic {
public:

    hydrogenic_11(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_11_x : public hydrogenic {
public:

    dell_hydrogenic_11_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_11_y : public hydrogenic {
public:

    dell_hydrogenic_11_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_11_z : public hydrogenic {
public:

    dell_hydrogenic_11_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_11 : public hydrogenic {
public:

    lapl_hydrogenic_11(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 11  -------------------------
*/

class hydrogenic_12 : public hydrogenic {
public:

    hydrogenic_12(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_12_x : public hydrogenic {
public:

    dell_hydrogenic_12_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_12_y : public hydrogenic {
public:

    dell_hydrogenic_12_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_12_z : public hydrogenic {
public:

    dell_hydrogenic_12_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_12 : public hydrogenic {
public:

    lapl_hydrogenic_12(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 12  -------------------------
*/

class hydrogenic_13 : public hydrogenic {
public:

    hydrogenic_13(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_13_x : public hydrogenic {
public:

    dell_hydrogenic_13_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_13_y : public hydrogenic {
public:

    dell_hydrogenic_13_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_13_z : public hydrogenic {
public:

    dell_hydrogenic_13_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_13 : public hydrogenic {
public:

    lapl_hydrogenic_13(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 13  -------------------------
*/

class hydrogenic_14 : public hydrogenic {
public:

    hydrogenic_14(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_14_x : public hydrogenic {
public:

    dell_hydrogenic_14_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_14_y : public hydrogenic {
public:

    dell_hydrogenic_14_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_14_z : public hydrogenic {
public:

    dell_hydrogenic_14_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_14 : public hydrogenic {
public:

    lapl_hydrogenic_14(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 14  -------------------------
*/

class hydrogenic_15 : public hydrogenic {
public:

    hydrogenic_15(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_15_x : public hydrogenic {
public:

    dell_hydrogenic_15_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_15_y : public hydrogenic {
public:

    dell_hydrogenic_15_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_15_z : public hydrogenic {
public:

    dell_hydrogenic_15_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_15 : public hydrogenic {
public:

    lapl_hydrogenic_15(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 15  -------------------------
*/

class hydrogenic_16 : public hydrogenic {
public:

    hydrogenic_16(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_16_x : public hydrogenic {
public:

    dell_hydrogenic_16_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_16_y : public hydrogenic {
public:

    dell_hydrogenic_16_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_16_z : public hydrogenic {
public:

    dell_hydrogenic_16_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_16 : public hydrogenic {
public:

    lapl_hydrogenic_16(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 16  -------------------------
*/

class hydrogenic_17 : public hydrogenic {
public:

    hydrogenic_17(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_17_x : public hydrogenic {
public:

    dell_hydrogenic_17_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_17_y : public hydrogenic {
public:

    dell_hydrogenic_17_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_hydrogenic_17_z : public hydrogenic {
public:

    dell_hydrogenic_17_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_hydrogenic_17 : public hydrogenic {
public:

    lapl_hydrogenic_17(double* k, double* k2, double* r22d, double* r2d, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 17  -------------------------
*/

