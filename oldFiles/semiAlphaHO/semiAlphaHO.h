/* 
 * File:   semiAlphaHO.h
 * Author: jorgehog
 *
 * Created on 27. juni 2012, 13:17
 */

#ifndef SEMIALPHAHO_H
#define	SEMIALPHAHO_H

class HO_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_1(double* alpha, double w);

    virtual double eval(const Walker* walker, int i);
};

class HO_2 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_2(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);
    
};

class HO_3 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_3(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_4 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_4(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_5 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_5(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_6 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_6(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_7 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_7(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_8 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_8(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_9 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_9(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);
  
};

class HO_10 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_10(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_11 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_11(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_12 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_12(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_13 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_13(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_14 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_14(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class HO_15 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    HO_15(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};


class dell_HO_1_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_1_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_1_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_1_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_2_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_2_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_2_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_2_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_2 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_2(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_3_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_3_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_3_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_3_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_3 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_3(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_4_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_4_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_4_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_4_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_4 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_4(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_5_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_5_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_5_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_5_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_5 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_5(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_6_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_6_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_6_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_6_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_6 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_6(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_7_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_7_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_7_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_7_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_7 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_7(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_8_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_8_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_8_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_8_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_8 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_8(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_9_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_9_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_9_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_9_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_9 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_9(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_10_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_10_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_10_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_10_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_10 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_10(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_11_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_11_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_11_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_11_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_11 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_11(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_12_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_12_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_12_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_12_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_12 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_12(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_13_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_13_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_13_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_13_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_13 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_13(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_14_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_14_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_14_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_14_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_14 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_14(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_15_0 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_15_0(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class dell_HO_15_1 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    dell_HO_15_1(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

class lapl_HO_15 : public BasisFunctions {
protected:
    double* alpha;
    double w;

public:

    lapl_HO_15(double* alpha, double w);
    
    virtual double eval(const Walker* walker, int i);

};

#endif	/* SEMIALPHAHO_H */