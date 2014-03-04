#pragma once

#include "../BasisFunctions.h"

//Superclass 
class HarmonicOscillator : public BasisFunctions {
protected:

    double* k;
    double* k2;
    double* exp_factor;

    double H;
    double x;
    double y;
    double x2;
    double y2;

public:

    HarmonicOscillator(double* k, double* k2, double* exp_factor);

};


/*
    Subclasses
*/


class HarmonicOscillator_0 : public HarmonicOscillator {
public:

    HarmonicOscillator_0(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_0_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_0_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_0_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_0_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_0 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_0(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 0  -------------------------
*/

class HarmonicOscillator_1 : public HarmonicOscillator {
public:

    HarmonicOscillator_1(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_1_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_1_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_1_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_1_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_1 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_1(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 1  -------------------------
*/

class HarmonicOscillator_2 : public HarmonicOscillator {
public:

    HarmonicOscillator_2(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_2_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_2_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_2_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_2_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_2 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_2(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 2  -------------------------
*/

class HarmonicOscillator_3 : public HarmonicOscillator {
public:

    HarmonicOscillator_3(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_3_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_3_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_3_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_3_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_3 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_3(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 3  -------------------------
*/

class HarmonicOscillator_4 : public HarmonicOscillator {
public:

    HarmonicOscillator_4(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_4_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_4_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_4_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_4_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_4 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_4(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 4  -------------------------
*/

class HarmonicOscillator_5 : public HarmonicOscillator {
public:

    HarmonicOscillator_5(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_5_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_5_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_5_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_5_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_5 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_5(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 5  -------------------------
*/

class HarmonicOscillator_6 : public HarmonicOscillator {
public:

    HarmonicOscillator_6(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_6_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_6_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_6_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_6_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_6 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_6(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 6  -------------------------
*/

class HarmonicOscillator_7 : public HarmonicOscillator {
public:

    HarmonicOscillator_7(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_7_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_7_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_7_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_7_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_7 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_7(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 7  -------------------------
*/

class HarmonicOscillator_8 : public HarmonicOscillator {
public:

    HarmonicOscillator_8(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_8_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_8_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_8_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_8_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_8 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_8(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 8  -------------------------
*/

class HarmonicOscillator_9 : public HarmonicOscillator {
public:

    HarmonicOscillator_9(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_9_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_9_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_9_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_9_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_9 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_9(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 9  -------------------------
*/

class HarmonicOscillator_10 : public HarmonicOscillator {
public:

    HarmonicOscillator_10(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_10_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_10_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_10_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_10_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_10 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_10(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 10  -------------------------
*/

class HarmonicOscillator_11 : public HarmonicOscillator {
public:

    HarmonicOscillator_11(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_11_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_11_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_11_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_11_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_11 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_11(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 11  -------------------------
*/

class HarmonicOscillator_12 : public HarmonicOscillator {
public:

    HarmonicOscillator_12(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_12_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_12_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_12_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_12_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_12 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_12(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 12  -------------------------
*/

class HarmonicOscillator_13 : public HarmonicOscillator {
public:

    HarmonicOscillator_13(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_13_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_13_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_13_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_13_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_13 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_13(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 13  -------------------------
*/

class HarmonicOscillator_14 : public HarmonicOscillator {
public:

    HarmonicOscillator_14(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_14_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_14_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_14_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_14_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_14 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_14(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 14  -------------------------
*/

class HarmonicOscillator_15 : public HarmonicOscillator {
public:

    HarmonicOscillator_15(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_15_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_15_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_15_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_15_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_15 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_15(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 15  -------------------------
*/

class HarmonicOscillator_16 : public HarmonicOscillator {
public:

    HarmonicOscillator_16(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_16_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_16_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_16_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_16_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_16 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_16(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 16  -------------------------
*/

class HarmonicOscillator_17 : public HarmonicOscillator {
public:

    HarmonicOscillator_17(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_17_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_17_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_17_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_17_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_17 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_17(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 17  -------------------------
*/

class HarmonicOscillator_18 : public HarmonicOscillator {
public:

    HarmonicOscillator_18(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_18_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_18_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_18_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_18_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_18 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_18(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 18  -------------------------
*/

class HarmonicOscillator_19 : public HarmonicOscillator {
public:

    HarmonicOscillator_19(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_19_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_19_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_19_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_19_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_19 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_19(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 19  -------------------------
*/

class HarmonicOscillator_20 : public HarmonicOscillator {
public:

    HarmonicOscillator_20(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_20_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_20_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_20_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_20_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_20 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_20(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 20  -------------------------
*/

class HarmonicOscillator_21 : public HarmonicOscillator {
public:

    HarmonicOscillator_21(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_21_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_21_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_21_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_21_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_21 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_21(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 21  -------------------------
*/

class HarmonicOscillator_22 : public HarmonicOscillator {
public:

    HarmonicOscillator_22(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_22_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_22_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_22_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_22_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_22 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_22(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 22  -------------------------
*/

class HarmonicOscillator_23 : public HarmonicOscillator {
public:

    HarmonicOscillator_23(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_23_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_23_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_23_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_23_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_23 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_23(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 23  -------------------------
*/

class HarmonicOscillator_24 : public HarmonicOscillator {
public:

    HarmonicOscillator_24(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_24_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_24_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_24_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_24_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_24 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_24(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 24  -------------------------
*/

class HarmonicOscillator_25 : public HarmonicOscillator {
public:

    HarmonicOscillator_25(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_25_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_25_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_25_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_25_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_25 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_25(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 25  -------------------------
*/

class HarmonicOscillator_26 : public HarmonicOscillator {
public:

    HarmonicOscillator_26(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_26_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_26_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_26_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_26_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_26 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_26(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 26  -------------------------
*/

class HarmonicOscillator_27 : public HarmonicOscillator {
public:

    HarmonicOscillator_27(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_27_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_27_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator_27_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_27_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator_27 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_27(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 27  -------------------------
*/

