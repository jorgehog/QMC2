
        
#ifndef HARMONICOSCILLATOR3D_H
#define HARMONICOSCILLATOR3D_H 

//Superclass 

class HarmonicOscillator3D : public BasisFunctions {
protected:

    double* k;
    double* k2;
    double* exp_factor;

    double H;
    double x;
    double y;
    double z;
    double x2;
    double y2;
    double z2;

public:

    HarmonicOscillator3D(double* k, double* k2, double* exp_factor);

};


/*
    Subclasses
*/


class HarmonicOscillator3D_0 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_0(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_0_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_0_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_0_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_0_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_0_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_0_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_0 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_0(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 0  -------------------------
*/

class HarmonicOscillator3D_1 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_1(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_1_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_1_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_1_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_1_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_1_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_1_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_1 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_1(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 1  -------------------------
*/

class HarmonicOscillator3D_2 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_2(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_2_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_2_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_2_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_2_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_2_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_2_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_2 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_2(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 2  -------------------------
*/

class HarmonicOscillator3D_3 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_3(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_3_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_3_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_3_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_3_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_3_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_3_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_3 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_3(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 3  -------------------------
*/

class HarmonicOscillator3D_4 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_4(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_4_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_4_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_4_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_4_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_4_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_4_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_4 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_4(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 4  -------------------------
*/

class HarmonicOscillator3D_5 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_5(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_5_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_5_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_5_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_5_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_5_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_5_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_5 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_5(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 5  -------------------------
*/

class HarmonicOscillator3D_6 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_6(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_6_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_6_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_6_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_6_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_6_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_6_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_6 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_6(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 6  -------------------------
*/

class HarmonicOscillator3D_7 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_7(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_7_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_7_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_7_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_7_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_7_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_7_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_7 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_7(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 7  -------------------------
*/

class HarmonicOscillator3D_8 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_8(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_8_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_8_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_8_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_8_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_8_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_8_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_8 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_8(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 8  -------------------------
*/

class HarmonicOscillator3D_9 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_9(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_9_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_9_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_9_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_9_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_9_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_9_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_9 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_9(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 9  -------------------------
*/

class HarmonicOscillator3D_10 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_10(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_10_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_10_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_10_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_10_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_10_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_10_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_10 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_10(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 10  -------------------------
*/

class HarmonicOscillator3D_11 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_11(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_11_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_11_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_11_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_11_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_11_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_11_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_11 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_11(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 11  -------------------------
*/

class HarmonicOscillator3D_12 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_12(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_12_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_12_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_12_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_12_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_12_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_12_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_12 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_12(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 12  -------------------------
*/

class HarmonicOscillator3D_13 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_13(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_13_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_13_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_13_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_13_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_13_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_13_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_13 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_13(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 13  -------------------------
*/

class HarmonicOscillator3D_14 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_14(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_14_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_14_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_14_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_14_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_14_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_14_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_14 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_14(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 14  -------------------------
*/

class HarmonicOscillator3D_15 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_15(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_15_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_15_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_15_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_15_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_15_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_15_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_15 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_15(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 15  -------------------------
*/

class HarmonicOscillator3D_16 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_16(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_16_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_16_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_16_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_16_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_16_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_16_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_16 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_16(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 16  -------------------------
*/

class HarmonicOscillator3D_17 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_17(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_17_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_17_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_17_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_17_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_17_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_17_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_17 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_17(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 17  -------------------------
*/

class HarmonicOscillator3D_18 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_18(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_18_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_18_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_18_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_18_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_18_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_18_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_18 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_18(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 18  -------------------------
*/

class HarmonicOscillator3D_19 : public HarmonicOscillator3D {
public:

    HarmonicOscillator3D_19(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_19_x : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_19_x(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_19_y : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_19_y(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class dell_HarmonicOscillator3D_19_z : public HarmonicOscillator3D {
public:

    dell_HarmonicOscillator3D_19_z(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

class lapl_HarmonicOscillator3D_19 : public HarmonicOscillator3D {
public:

    lapl_HarmonicOscillator3D_19(double* k, double* k2, double* exp_factor);
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 19  -------------------------
*/


#endif /* HARMONICOSCILLATOR3D_H */
