#ifndef GAUSSIANS_H
#define GAUSSIANS_H 

#include "../BasisFunctions.h"

//Superclass 

class gaussians : public BasisFunctions {
protected:

    double* exp_factor;

    double a;
    double P;
    double x;
    double y;
    double z;
    double x2;
    double y2;
    double z2;

public:

    gaussians(double* exp_factor, double a);

};


/*
    Subclasses
*/


class gaussians_000 : public gaussians {
public:

    gaussians_000(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_000_x : public gaussians {
public:

    dell_gaussians_000_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_000_y : public gaussians {
public:

    dell_gaussians_000_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_000_z : public gaussians {
public:

    dell_gaussians_000_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_000 : public gaussians {
public:

    lapl_gaussians_000(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 0  -------------------------
*/

class gaussians_001 : public gaussians {
public:

    gaussians_001(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_001_x : public gaussians {
public:

    dell_gaussians_001_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_001_y : public gaussians {
public:

    dell_gaussians_001_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_001_z : public gaussians {
public:

    dell_gaussians_001_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_001 : public gaussians {
public:

    lapl_gaussians_001(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 1  -------------------------
*/

class gaussians_010 : public gaussians {
public:

    gaussians_010(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_010_x : public gaussians {
public:

    dell_gaussians_010_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_010_y : public gaussians {
public:

    dell_gaussians_010_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_010_z : public gaussians {
public:

    dell_gaussians_010_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_010 : public gaussians {
public:

    lapl_gaussians_010(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 2  -------------------------
*/

class gaussians_100 : public gaussians {
public:

    gaussians_100(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_100_x : public gaussians {
public:

    dell_gaussians_100_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_100_y : public gaussians {
public:

    dell_gaussians_100_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_100_z : public gaussians {
public:

    dell_gaussians_100_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_100 : public gaussians {
public:

    lapl_gaussians_100(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 3  -------------------------
*/

class gaussians_002 : public gaussians {
public:

    gaussians_002(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_002_x : public gaussians {
public:

    dell_gaussians_002_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_002_y : public gaussians {
public:

    dell_gaussians_002_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_002_z : public gaussians {
public:

    dell_gaussians_002_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_002 : public gaussians {
public:

    lapl_gaussians_002(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 4  -------------------------
*/

class gaussians_011 : public gaussians {
public:

    gaussians_011(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_011_x : public gaussians {
public:

    dell_gaussians_011_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_011_y : public gaussians {
public:

    dell_gaussians_011_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_011_z : public gaussians {
public:

    dell_gaussians_011_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_011 : public gaussians {
public:

    lapl_gaussians_011(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 5  -------------------------
*/

class gaussians_020 : public gaussians {
public:

    gaussians_020(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_020_x : public gaussians {
public:

    dell_gaussians_020_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_020_y : public gaussians {
public:

    dell_gaussians_020_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_020_z : public gaussians {
public:

    dell_gaussians_020_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_020 : public gaussians {
public:

    lapl_gaussians_020(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 6  -------------------------
*/

class gaussians_101 : public gaussians {
public:

    gaussians_101(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_101_x : public gaussians {
public:

    dell_gaussians_101_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_101_y : public gaussians {
public:

    dell_gaussians_101_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_101_z : public gaussians {
public:

    dell_gaussians_101_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_101 : public gaussians {
public:

    lapl_gaussians_101(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 7  -------------------------
*/

class gaussians_110 : public gaussians {
public:

    gaussians_110(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_110_x : public gaussians {
public:

    dell_gaussians_110_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_110_y : public gaussians {
public:

    dell_gaussians_110_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_110_z : public gaussians {
public:

    dell_gaussians_110_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_110 : public gaussians {
public:

    lapl_gaussians_110(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 8  -------------------------
*/

class gaussians_200 : public gaussians {
public:

    gaussians_200(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_200_x : public gaussians {
public:

    dell_gaussians_200_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_200_y : public gaussians {
public:

    dell_gaussians_200_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_200_z : public gaussians {
public:

    dell_gaussians_200_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_200 : public gaussians {
public:

    lapl_gaussians_200(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 9  -------------------------
*/

class gaussians_003 : public gaussians {
public:

    gaussians_003(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_003_x : public gaussians {
public:

    dell_gaussians_003_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_003_y : public gaussians {
public:

    dell_gaussians_003_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_003_z : public gaussians {
public:

    dell_gaussians_003_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_003 : public gaussians {
public:

    lapl_gaussians_003(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 10  -------------------------
*/

class gaussians_012 : public gaussians {
public:

    gaussians_012(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_012_x : public gaussians {
public:

    dell_gaussians_012_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_012_y : public gaussians {
public:

    dell_gaussians_012_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_012_z : public gaussians {
public:

    dell_gaussians_012_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_012 : public gaussians {
public:

    lapl_gaussians_012(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 11  -------------------------
*/

class gaussians_021 : public gaussians {
public:

    gaussians_021(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_021_x : public gaussians {
public:

    dell_gaussians_021_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_021_y : public gaussians {
public:

    dell_gaussians_021_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_021_z : public gaussians {
public:

    dell_gaussians_021_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_021 : public gaussians {
public:

    lapl_gaussians_021(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 12  -------------------------
*/

class gaussians_030 : public gaussians {
public:

    gaussians_030(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_030_x : public gaussians {
public:

    dell_gaussians_030_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_030_y : public gaussians {
public:

    dell_gaussians_030_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_030_z : public gaussians {
public:

    dell_gaussians_030_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_030 : public gaussians {
public:

    lapl_gaussians_030(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 13  -------------------------
*/

class gaussians_102 : public gaussians {
public:

    gaussians_102(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_102_x : public gaussians {
public:

    dell_gaussians_102_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_102_y : public gaussians {
public:

    dell_gaussians_102_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_102_z : public gaussians {
public:

    dell_gaussians_102_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_102 : public gaussians {
public:

    lapl_gaussians_102(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 14  -------------------------
*/

class gaussians_111 : public gaussians {
public:

    gaussians_111(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_111_x : public gaussians {
public:

    dell_gaussians_111_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_111_y : public gaussians {
public:

    dell_gaussians_111_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_111_z : public gaussians {
public:

    dell_gaussians_111_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_111 : public gaussians {
public:

    lapl_gaussians_111(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 15  -------------------------
*/

class gaussians_120 : public gaussians {
public:

    gaussians_120(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_120_x : public gaussians {
public:

    dell_gaussians_120_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_120_y : public gaussians {
public:

    dell_gaussians_120_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_120_z : public gaussians {
public:

    dell_gaussians_120_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_120 : public gaussians {
public:

    lapl_gaussians_120(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 16  -------------------------
*/

class gaussians_201 : public gaussians {
public:

    gaussians_201(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_201_x : public gaussians {
public:

    dell_gaussians_201_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_201_y : public gaussians {
public:

    dell_gaussians_201_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_201_z : public gaussians {
public:

    dell_gaussians_201_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_201 : public gaussians {
public:

    lapl_gaussians_201(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 17  -------------------------
*/

class gaussians_210 : public gaussians {
public:

    gaussians_210(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_210_x : public gaussians {
public:

    dell_gaussians_210_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_210_y : public gaussians {
public:

    dell_gaussians_210_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_210_z : public gaussians {
public:

    dell_gaussians_210_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_210 : public gaussians {
public:

    lapl_gaussians_210(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 18  -------------------------
*/

class gaussians_300 : public gaussians {
public:

    gaussians_300(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_300_x : public gaussians {
public:

    dell_gaussians_300_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_300_y : public gaussians {
public:

    dell_gaussians_300_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_300_z : public gaussians {
public:

    dell_gaussians_300_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_300 : public gaussians {
public:

    lapl_gaussians_300(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 19  -------------------------
*/

class gaussians_004 : public gaussians {
public:

    gaussians_004(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_004_x : public gaussians {
public:

    dell_gaussians_004_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_004_y : public gaussians {
public:

    dell_gaussians_004_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_004_z : public gaussians {
public:

    dell_gaussians_004_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_004 : public gaussians {
public:

    lapl_gaussians_004(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 20  -------------------------
*/

class gaussians_013 : public gaussians {
public:

    gaussians_013(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_013_x : public gaussians {
public:

    dell_gaussians_013_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_013_y : public gaussians {
public:

    dell_gaussians_013_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_013_z : public gaussians {
public:

    dell_gaussians_013_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_013 : public gaussians {
public:

    lapl_gaussians_013(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 21  -------------------------
*/

class gaussians_022 : public gaussians {
public:

    gaussians_022(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_022_x : public gaussians {
public:

    dell_gaussians_022_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_022_y : public gaussians {
public:

    dell_gaussians_022_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_022_z : public gaussians {
public:

    dell_gaussians_022_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_022 : public gaussians {
public:

    lapl_gaussians_022(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 22  -------------------------
*/

class gaussians_031 : public gaussians {
public:

    gaussians_031(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_031_x : public gaussians {
public:

    dell_gaussians_031_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_031_y : public gaussians {
public:

    dell_gaussians_031_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_031_z : public gaussians {
public:

    dell_gaussians_031_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_031 : public gaussians {
public:

    lapl_gaussians_031(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 23  -------------------------
*/

class gaussians_040 : public gaussians {
public:

    gaussians_040(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_040_x : public gaussians {
public:

    dell_gaussians_040_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_040_y : public gaussians {
public:

    dell_gaussians_040_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_040_z : public gaussians {
public:

    dell_gaussians_040_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_040 : public gaussians {
public:

    lapl_gaussians_040(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 24  -------------------------
*/

class gaussians_103 : public gaussians {
public:

    gaussians_103(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_103_x : public gaussians {
public:

    dell_gaussians_103_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_103_y : public gaussians {
public:

    dell_gaussians_103_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_103_z : public gaussians {
public:

    dell_gaussians_103_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_103 : public gaussians {
public:

    lapl_gaussians_103(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 25  -------------------------
*/

class gaussians_112 : public gaussians {
public:

    gaussians_112(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_112_x : public gaussians {
public:

    dell_gaussians_112_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_112_y : public gaussians {
public:

    dell_gaussians_112_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_112_z : public gaussians {
public:

    dell_gaussians_112_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_112 : public gaussians {
public:

    lapl_gaussians_112(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 26  -------------------------
*/

class gaussians_121 : public gaussians {
public:

    gaussians_121(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_121_x : public gaussians {
public:

    dell_gaussians_121_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_121_y : public gaussians {
public:

    dell_gaussians_121_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_121_z : public gaussians {
public:

    dell_gaussians_121_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_121 : public gaussians {
public:

    lapl_gaussians_121(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 27  -------------------------
*/

class gaussians_130 : public gaussians {
public:

    gaussians_130(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_130_x : public gaussians {
public:

    dell_gaussians_130_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_130_y : public gaussians {
public:

    dell_gaussians_130_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_130_z : public gaussians {
public:

    dell_gaussians_130_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_130 : public gaussians {
public:

    lapl_gaussians_130(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 28  -------------------------
*/

class gaussians_202 : public gaussians {
public:

    gaussians_202(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_202_x : public gaussians {
public:

    dell_gaussians_202_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_202_y : public gaussians {
public:

    dell_gaussians_202_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_202_z : public gaussians {
public:

    dell_gaussians_202_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_202 : public gaussians {
public:

    lapl_gaussians_202(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 29  -------------------------
*/

class gaussians_211 : public gaussians {
public:

    gaussians_211(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_211_x : public gaussians {
public:

    dell_gaussians_211_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_211_y : public gaussians {
public:

    dell_gaussians_211_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_211_z : public gaussians {
public:

    dell_gaussians_211_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_211 : public gaussians {
public:

    lapl_gaussians_211(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 30  -------------------------
*/

class gaussians_220 : public gaussians {
public:

    gaussians_220(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_220_x : public gaussians {
public:

    dell_gaussians_220_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_220_y : public gaussians {
public:

    dell_gaussians_220_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_220_z : public gaussians {
public:

    dell_gaussians_220_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_220 : public gaussians {
public:

    lapl_gaussians_220(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 31  -------------------------
*/

class gaussians_301 : public gaussians {
public:

    gaussians_301(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_301_x : public gaussians {
public:

    dell_gaussians_301_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_301_y : public gaussians {
public:

    dell_gaussians_301_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_301_z : public gaussians {
public:

    dell_gaussians_301_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_301 : public gaussians {
public:

    lapl_gaussians_301(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 32  -------------------------
*/

class gaussians_310 : public gaussians {
public:

    gaussians_310(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_310_x : public gaussians {
public:

    dell_gaussians_310_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_310_y : public gaussians {
public:

    dell_gaussians_310_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_310_z : public gaussians {
public:

    dell_gaussians_310_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_310 : public gaussians {
public:

    lapl_gaussians_310(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 33  -------------------------
*/

class gaussians_400 : public gaussians {
public:

    gaussians_400(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_400_x : public gaussians {
public:

    dell_gaussians_400_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_400_y : public gaussians {
public:

    dell_gaussians_400_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_400_z : public gaussians {
public:

    dell_gaussians_400_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_400 : public gaussians {
public:

    lapl_gaussians_400(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 34  -------------------------
*/

class gaussians_005 : public gaussians {
public:

    gaussians_005(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_005_x : public gaussians {
public:

    dell_gaussians_005_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_005_y : public gaussians {
public:

    dell_gaussians_005_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_005_z : public gaussians {
public:

    dell_gaussians_005_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_005 : public gaussians {
public:

    lapl_gaussians_005(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 35  -------------------------
*/

class gaussians_014 : public gaussians {
public:

    gaussians_014(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_014_x : public gaussians {
public:

    dell_gaussians_014_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_014_y : public gaussians {
public:

    dell_gaussians_014_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_014_z : public gaussians {
public:

    dell_gaussians_014_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_014 : public gaussians {
public:

    lapl_gaussians_014(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 36  -------------------------
*/

class gaussians_023 : public gaussians {
public:

    gaussians_023(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_023_x : public gaussians {
public:

    dell_gaussians_023_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_023_y : public gaussians {
public:

    dell_gaussians_023_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_023_z : public gaussians {
public:

    dell_gaussians_023_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_023 : public gaussians {
public:

    lapl_gaussians_023(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 37  -------------------------
*/

class gaussians_032 : public gaussians {
public:

    gaussians_032(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_032_x : public gaussians {
public:

    dell_gaussians_032_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_032_y : public gaussians {
public:

    dell_gaussians_032_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_032_z : public gaussians {
public:

    dell_gaussians_032_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_032 : public gaussians {
public:

    lapl_gaussians_032(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 38  -------------------------
*/

class gaussians_041 : public gaussians {
public:

    gaussians_041(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_041_x : public gaussians {
public:

    dell_gaussians_041_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_041_y : public gaussians {
public:

    dell_gaussians_041_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_041_z : public gaussians {
public:

    dell_gaussians_041_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_041 : public gaussians {
public:

    lapl_gaussians_041(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 39  -------------------------
*/

class gaussians_050 : public gaussians {
public:

    gaussians_050(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_050_x : public gaussians {
public:

    dell_gaussians_050_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_050_y : public gaussians {
public:

    dell_gaussians_050_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_050_z : public gaussians {
public:

    dell_gaussians_050_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_050 : public gaussians {
public:

    lapl_gaussians_050(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 40  -------------------------
*/

class gaussians_104 : public gaussians {
public:

    gaussians_104(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_104_x : public gaussians {
public:

    dell_gaussians_104_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_104_y : public gaussians {
public:

    dell_gaussians_104_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_104_z : public gaussians {
public:

    dell_gaussians_104_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_104 : public gaussians {
public:

    lapl_gaussians_104(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 41  -------------------------
*/

class gaussians_113 : public gaussians {
public:

    gaussians_113(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_113_x : public gaussians {
public:

    dell_gaussians_113_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_113_y : public gaussians {
public:

    dell_gaussians_113_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_113_z : public gaussians {
public:

    dell_gaussians_113_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_113 : public gaussians {
public:

    lapl_gaussians_113(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 42  -------------------------
*/

class gaussians_122 : public gaussians {
public:

    gaussians_122(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_122_x : public gaussians {
public:

    dell_gaussians_122_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_122_y : public gaussians {
public:

    dell_gaussians_122_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_122_z : public gaussians {
public:

    dell_gaussians_122_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_122 : public gaussians {
public:

    lapl_gaussians_122(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 43  -------------------------
*/

class gaussians_131 : public gaussians {
public:

    gaussians_131(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_131_x : public gaussians {
public:

    dell_gaussians_131_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_131_y : public gaussians {
public:

    dell_gaussians_131_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_131_z : public gaussians {
public:

    dell_gaussians_131_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_131 : public gaussians {
public:

    lapl_gaussians_131(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 44  -------------------------
*/

class gaussians_140 : public gaussians {
public:

    gaussians_140(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_140_x : public gaussians {
public:

    dell_gaussians_140_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_140_y : public gaussians {
public:

    dell_gaussians_140_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_140_z : public gaussians {
public:

    dell_gaussians_140_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_140 : public gaussians {
public:

    lapl_gaussians_140(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 45  -------------------------
*/

class gaussians_203 : public gaussians {
public:

    gaussians_203(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_203_x : public gaussians {
public:

    dell_gaussians_203_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_203_y : public gaussians {
public:

    dell_gaussians_203_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_203_z : public gaussians {
public:

    dell_gaussians_203_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_203 : public gaussians {
public:

    lapl_gaussians_203(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 46  -------------------------
*/

class gaussians_212 : public gaussians {
public:

    gaussians_212(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_212_x : public gaussians {
public:

    dell_gaussians_212_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_212_y : public gaussians {
public:

    dell_gaussians_212_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_212_z : public gaussians {
public:

    dell_gaussians_212_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_212 : public gaussians {
public:

    lapl_gaussians_212(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 47  -------------------------
*/

class gaussians_221 : public gaussians {
public:

    gaussians_221(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_221_x : public gaussians {
public:

    dell_gaussians_221_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_221_y : public gaussians {
public:

    dell_gaussians_221_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_221_z : public gaussians {
public:

    dell_gaussians_221_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_221 : public gaussians {
public:

    lapl_gaussians_221(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 48  -------------------------
*/

class gaussians_230 : public gaussians {
public:

    gaussians_230(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_230_x : public gaussians {
public:

    dell_gaussians_230_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_230_y : public gaussians {
public:

    dell_gaussians_230_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_230_z : public gaussians {
public:

    dell_gaussians_230_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_230 : public gaussians {
public:

    lapl_gaussians_230(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 49  -------------------------
*/

class gaussians_302 : public gaussians {
public:

    gaussians_302(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_302_x : public gaussians {
public:

    dell_gaussians_302_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_302_y : public gaussians {
public:

    dell_gaussians_302_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_302_z : public gaussians {
public:

    dell_gaussians_302_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_302 : public gaussians {
public:

    lapl_gaussians_302(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 50  -------------------------
*/

class gaussians_311 : public gaussians {
public:

    gaussians_311(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_311_x : public gaussians {
public:

    dell_gaussians_311_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_311_y : public gaussians {
public:

    dell_gaussians_311_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_311_z : public gaussians {
public:

    dell_gaussians_311_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_311 : public gaussians {
public:

    lapl_gaussians_311(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 51  -------------------------
*/

class gaussians_320 : public gaussians {
public:

    gaussians_320(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_320_x : public gaussians {
public:

    dell_gaussians_320_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_320_y : public gaussians {
public:

    dell_gaussians_320_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_320_z : public gaussians {
public:

    dell_gaussians_320_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_320 : public gaussians {
public:

    lapl_gaussians_320(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 52  -------------------------
*/

class gaussians_401 : public gaussians {
public:

    gaussians_401(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_401_x : public gaussians {
public:

    dell_gaussians_401_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_401_y : public gaussians {
public:

    dell_gaussians_401_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_401_z : public gaussians {
public:

    dell_gaussians_401_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_401 : public gaussians {
public:

    lapl_gaussians_401(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 53  -------------------------
*/

class gaussians_410 : public gaussians {
public:

    gaussians_410(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_410_x : public gaussians {
public:

    dell_gaussians_410_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_410_y : public gaussians {
public:

    dell_gaussians_410_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_410_z : public gaussians {
public:

    dell_gaussians_410_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_410 : public gaussians {
public:

    lapl_gaussians_410(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 54  -------------------------
*/

class gaussians_500 : public gaussians {
public:

    gaussians_500(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_500_x : public gaussians {
public:

    dell_gaussians_500_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_500_y : public gaussians {
public:

    dell_gaussians_500_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_500_z : public gaussians {
public:

    dell_gaussians_500_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_500 : public gaussians {
public:

    lapl_gaussians_500(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 55  -------------------------
*/

class gaussians_006 : public gaussians {
public:

    gaussians_006(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_006_x : public gaussians {
public:

    dell_gaussians_006_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_006_y : public gaussians {
public:

    dell_gaussians_006_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_006_z : public gaussians {
public:

    dell_gaussians_006_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_006 : public gaussians {
public:

    lapl_gaussians_006(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 56  -------------------------
*/

class gaussians_015 : public gaussians {
public:

    gaussians_015(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_015_x : public gaussians {
public:

    dell_gaussians_015_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_015_y : public gaussians {
public:

    dell_gaussians_015_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_015_z : public gaussians {
public:

    dell_gaussians_015_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_015 : public gaussians {
public:

    lapl_gaussians_015(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 57  -------------------------
*/

class gaussians_024 : public gaussians {
public:

    gaussians_024(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_024_x : public gaussians {
public:

    dell_gaussians_024_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_024_y : public gaussians {
public:

    dell_gaussians_024_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_024_z : public gaussians {
public:

    dell_gaussians_024_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_024 : public gaussians {
public:

    lapl_gaussians_024(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 58  -------------------------
*/

class gaussians_033 : public gaussians {
public:

    gaussians_033(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_033_x : public gaussians {
public:

    dell_gaussians_033_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_033_y : public gaussians {
public:

    dell_gaussians_033_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_033_z : public gaussians {
public:

    dell_gaussians_033_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_033 : public gaussians {
public:

    lapl_gaussians_033(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 59  -------------------------
*/

class gaussians_042 : public gaussians {
public:

    gaussians_042(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_042_x : public gaussians {
public:

    dell_gaussians_042_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_042_y : public gaussians {
public:

    dell_gaussians_042_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_042_z : public gaussians {
public:

    dell_gaussians_042_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_042 : public gaussians {
public:

    lapl_gaussians_042(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 60  -------------------------
*/

class gaussians_051 : public gaussians {
public:

    gaussians_051(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_051_x : public gaussians {
public:

    dell_gaussians_051_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_051_y : public gaussians {
public:

    dell_gaussians_051_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_051_z : public gaussians {
public:

    dell_gaussians_051_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_051 : public gaussians {
public:

    lapl_gaussians_051(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 61  -------------------------
*/

class gaussians_060 : public gaussians {
public:

    gaussians_060(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_060_x : public gaussians {
public:

    dell_gaussians_060_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_060_y : public gaussians {
public:

    dell_gaussians_060_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_060_z : public gaussians {
public:

    dell_gaussians_060_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_060 : public gaussians {
public:

    lapl_gaussians_060(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 62  -------------------------
*/

class gaussians_105 : public gaussians {
public:

    gaussians_105(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_105_x : public gaussians {
public:

    dell_gaussians_105_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_105_y : public gaussians {
public:

    dell_gaussians_105_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_105_z : public gaussians {
public:

    dell_gaussians_105_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_105 : public gaussians {
public:

    lapl_gaussians_105(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 63  -------------------------
*/

class gaussians_114 : public gaussians {
public:

    gaussians_114(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_114_x : public gaussians {
public:

    dell_gaussians_114_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_114_y : public gaussians {
public:

    dell_gaussians_114_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_114_z : public gaussians {
public:

    dell_gaussians_114_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_114 : public gaussians {
public:

    lapl_gaussians_114(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 64  -------------------------
*/

class gaussians_123 : public gaussians {
public:

    gaussians_123(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_123_x : public gaussians {
public:

    dell_gaussians_123_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_123_y : public gaussians {
public:

    dell_gaussians_123_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_123_z : public gaussians {
public:

    dell_gaussians_123_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_123 : public gaussians {
public:

    lapl_gaussians_123(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 65  -------------------------
*/

class gaussians_132 : public gaussians {
public:

    gaussians_132(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_132_x : public gaussians {
public:

    dell_gaussians_132_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_132_y : public gaussians {
public:

    dell_gaussians_132_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_132_z : public gaussians {
public:

    dell_gaussians_132_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_132 : public gaussians {
public:

    lapl_gaussians_132(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 66  -------------------------
*/

class gaussians_141 : public gaussians {
public:

    gaussians_141(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_141_x : public gaussians {
public:

    dell_gaussians_141_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_141_y : public gaussians {
public:

    dell_gaussians_141_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_141_z : public gaussians {
public:

    dell_gaussians_141_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_141 : public gaussians {
public:

    lapl_gaussians_141(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 67  -------------------------
*/

class gaussians_150 : public gaussians {
public:

    gaussians_150(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_150_x : public gaussians {
public:

    dell_gaussians_150_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_150_y : public gaussians {
public:

    dell_gaussians_150_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_150_z : public gaussians {
public:

    dell_gaussians_150_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_150 : public gaussians {
public:

    lapl_gaussians_150(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 68  -------------------------
*/

class gaussians_204 : public gaussians {
public:

    gaussians_204(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_204_x : public gaussians {
public:

    dell_gaussians_204_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_204_y : public gaussians {
public:

    dell_gaussians_204_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_204_z : public gaussians {
public:

    dell_gaussians_204_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_204 : public gaussians {
public:

    lapl_gaussians_204(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 69  -------------------------
*/

class gaussians_213 : public gaussians {
public:

    gaussians_213(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_213_x : public gaussians {
public:

    dell_gaussians_213_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_213_y : public gaussians {
public:

    dell_gaussians_213_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_213_z : public gaussians {
public:

    dell_gaussians_213_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_213 : public gaussians {
public:

    lapl_gaussians_213(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 70  -------------------------
*/

class gaussians_222 : public gaussians {
public:

    gaussians_222(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_222_x : public gaussians {
public:

    dell_gaussians_222_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_222_y : public gaussians {
public:

    dell_gaussians_222_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_222_z : public gaussians {
public:

    dell_gaussians_222_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_222 : public gaussians {
public:

    lapl_gaussians_222(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 71  -------------------------
*/

class gaussians_231 : public gaussians {
public:

    gaussians_231(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_231_x : public gaussians {
public:

    dell_gaussians_231_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_231_y : public gaussians {
public:

    dell_gaussians_231_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_231_z : public gaussians {
public:

    dell_gaussians_231_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_231 : public gaussians {
public:

    lapl_gaussians_231(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 72  -------------------------
*/

class gaussians_240 : public gaussians {
public:

    gaussians_240(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_240_x : public gaussians {
public:

    dell_gaussians_240_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_240_y : public gaussians {
public:

    dell_gaussians_240_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_240_z : public gaussians {
public:

    dell_gaussians_240_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_240 : public gaussians {
public:

    lapl_gaussians_240(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 73  -------------------------
*/

class gaussians_303 : public gaussians {
public:

    gaussians_303(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_303_x : public gaussians {
public:

    dell_gaussians_303_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_303_y : public gaussians {
public:

    dell_gaussians_303_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_303_z : public gaussians {
public:

    dell_gaussians_303_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_303 : public gaussians {
public:

    lapl_gaussians_303(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 74  -------------------------
*/

class gaussians_312 : public gaussians {
public:

    gaussians_312(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_312_x : public gaussians {
public:

    dell_gaussians_312_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_312_y : public gaussians {
public:

    dell_gaussians_312_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_312_z : public gaussians {
public:

    dell_gaussians_312_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_312 : public gaussians {
public:

    lapl_gaussians_312(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 75  -------------------------
*/

class gaussians_321 : public gaussians {
public:

    gaussians_321(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_321_x : public gaussians {
public:

    dell_gaussians_321_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_321_y : public gaussians {
public:

    dell_gaussians_321_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_321_z : public gaussians {
public:

    dell_gaussians_321_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_321 : public gaussians {
public:

    lapl_gaussians_321(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 76  -------------------------
*/

class gaussians_330 : public gaussians {
public:

    gaussians_330(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_330_x : public gaussians {
public:

    dell_gaussians_330_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_330_y : public gaussians {
public:

    dell_gaussians_330_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_330_z : public gaussians {
public:

    dell_gaussians_330_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_330 : public gaussians {
public:

    lapl_gaussians_330(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 77  -------------------------
*/

class gaussians_402 : public gaussians {
public:

    gaussians_402(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_402_x : public gaussians {
public:

    dell_gaussians_402_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_402_y : public gaussians {
public:

    dell_gaussians_402_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_402_z : public gaussians {
public:

    dell_gaussians_402_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_402 : public gaussians {
public:

    lapl_gaussians_402(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 78  -------------------------
*/

class gaussians_411 : public gaussians {
public:

    gaussians_411(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_411_x : public gaussians {
public:

    dell_gaussians_411_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_411_y : public gaussians {
public:

    dell_gaussians_411_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_411_z : public gaussians {
public:

    dell_gaussians_411_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_411 : public gaussians {
public:

    lapl_gaussians_411(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 79  -------------------------
*/

class gaussians_420 : public gaussians {
public:

    gaussians_420(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_420_x : public gaussians {
public:

    dell_gaussians_420_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_420_y : public gaussians {
public:

    dell_gaussians_420_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_420_z : public gaussians {
public:

    dell_gaussians_420_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_420 : public gaussians {
public:

    lapl_gaussians_420(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 80  -------------------------
*/

class gaussians_501 : public gaussians {
public:

    gaussians_501(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_501_x : public gaussians {
public:

    dell_gaussians_501_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_501_y : public gaussians {
public:

    dell_gaussians_501_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_501_z : public gaussians {
public:

    dell_gaussians_501_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_501 : public gaussians {
public:

    lapl_gaussians_501(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 81  -------------------------
*/

class gaussians_510 : public gaussians {
public:

    gaussians_510(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_510_x : public gaussians {
public:

    dell_gaussians_510_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_510_y : public gaussians {
public:

    dell_gaussians_510_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_510_z : public gaussians {
public:

    dell_gaussians_510_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_510 : public gaussians {
public:

    lapl_gaussians_510(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 82  -------------------------
*/

class gaussians_600 : public gaussians {
public:

    gaussians_600(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_600_x : public gaussians {
public:

    dell_gaussians_600_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_600_y : public gaussians {
public:

    dell_gaussians_600_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_600_z : public gaussians {
public:

    dell_gaussians_600_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_600 : public gaussians {
public:

    lapl_gaussians_600(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 83  -------------------------
*/

class gaussians_007 : public gaussians {
public:

    gaussians_007(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_007_x : public gaussians {
public:

    dell_gaussians_007_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_007_y : public gaussians {
public:

    dell_gaussians_007_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_007_z : public gaussians {
public:

    dell_gaussians_007_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_007 : public gaussians {
public:

    lapl_gaussians_007(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 84  -------------------------
*/

class gaussians_016 : public gaussians {
public:

    gaussians_016(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_016_x : public gaussians {
public:

    dell_gaussians_016_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_016_y : public gaussians {
public:

    dell_gaussians_016_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_016_z : public gaussians {
public:

    dell_gaussians_016_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_016 : public gaussians {
public:

    lapl_gaussians_016(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 85  -------------------------
*/

class gaussians_025 : public gaussians {
public:

    gaussians_025(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_025_x : public gaussians {
public:

    dell_gaussians_025_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_025_y : public gaussians {
public:

    dell_gaussians_025_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_025_z : public gaussians {
public:

    dell_gaussians_025_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_025 : public gaussians {
public:

    lapl_gaussians_025(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 86  -------------------------
*/

class gaussians_034 : public gaussians {
public:

    gaussians_034(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_034_x : public gaussians {
public:

    dell_gaussians_034_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_034_y : public gaussians {
public:

    dell_gaussians_034_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_034_z : public gaussians {
public:

    dell_gaussians_034_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_034 : public gaussians {
public:

    lapl_gaussians_034(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 87  -------------------------
*/

class gaussians_043 : public gaussians {
public:

    gaussians_043(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_043_x : public gaussians {
public:

    dell_gaussians_043_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_043_y : public gaussians {
public:

    dell_gaussians_043_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_043_z : public gaussians {
public:

    dell_gaussians_043_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_043 : public gaussians {
public:

    lapl_gaussians_043(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 88  -------------------------
*/

class gaussians_052 : public gaussians {
public:

    gaussians_052(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_052_x : public gaussians {
public:

    dell_gaussians_052_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_052_y : public gaussians {
public:

    dell_gaussians_052_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_052_z : public gaussians {
public:

    dell_gaussians_052_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_052 : public gaussians {
public:

    lapl_gaussians_052(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 89  -------------------------
*/

class gaussians_061 : public gaussians {
public:

    gaussians_061(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_061_x : public gaussians {
public:

    dell_gaussians_061_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_061_y : public gaussians {
public:

    dell_gaussians_061_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_061_z : public gaussians {
public:

    dell_gaussians_061_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_061 : public gaussians {
public:

    lapl_gaussians_061(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 90  -------------------------
*/

class gaussians_070 : public gaussians {
public:

    gaussians_070(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_070_x : public gaussians {
public:

    dell_gaussians_070_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_070_y : public gaussians {
public:

    dell_gaussians_070_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_070_z : public gaussians {
public:

    dell_gaussians_070_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_070 : public gaussians {
public:

    lapl_gaussians_070(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 91  -------------------------
*/

class gaussians_106 : public gaussians {
public:

    gaussians_106(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_106_x : public gaussians {
public:

    dell_gaussians_106_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_106_y : public gaussians {
public:

    dell_gaussians_106_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_106_z : public gaussians {
public:

    dell_gaussians_106_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_106 : public gaussians {
public:

    lapl_gaussians_106(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 92  -------------------------
*/

class gaussians_115 : public gaussians {
public:

    gaussians_115(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_115_x : public gaussians {
public:

    dell_gaussians_115_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_115_y : public gaussians {
public:

    dell_gaussians_115_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_115_z : public gaussians {
public:

    dell_gaussians_115_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_115 : public gaussians {
public:

    lapl_gaussians_115(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 93  -------------------------
*/

class gaussians_124 : public gaussians {
public:

    gaussians_124(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_124_x : public gaussians {
public:

    dell_gaussians_124_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_124_y : public gaussians {
public:

    dell_gaussians_124_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_124_z : public gaussians {
public:

    dell_gaussians_124_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_124 : public gaussians {
public:

    lapl_gaussians_124(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 94  -------------------------
*/

class gaussians_133 : public gaussians {
public:

    gaussians_133(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_133_x : public gaussians {
public:

    dell_gaussians_133_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_133_y : public gaussians {
public:

    dell_gaussians_133_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_133_z : public gaussians {
public:

    dell_gaussians_133_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_133 : public gaussians {
public:

    lapl_gaussians_133(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 95  -------------------------
*/

class gaussians_142 : public gaussians {
public:

    gaussians_142(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_142_x : public gaussians {
public:

    dell_gaussians_142_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_142_y : public gaussians {
public:

    dell_gaussians_142_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_142_z : public gaussians {
public:

    dell_gaussians_142_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_142 : public gaussians {
public:

    lapl_gaussians_142(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 96  -------------------------
*/

class gaussians_151 : public gaussians {
public:

    gaussians_151(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_151_x : public gaussians {
public:

    dell_gaussians_151_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_151_y : public gaussians {
public:

    dell_gaussians_151_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_151_z : public gaussians {
public:

    dell_gaussians_151_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_151 : public gaussians {
public:

    lapl_gaussians_151(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 97  -------------------------
*/

class gaussians_160 : public gaussians {
public:

    gaussians_160(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_160_x : public gaussians {
public:

    dell_gaussians_160_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_160_y : public gaussians {
public:

    dell_gaussians_160_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_160_z : public gaussians {
public:

    dell_gaussians_160_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_160 : public gaussians {
public:

    lapl_gaussians_160(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 98  -------------------------
*/

class gaussians_205 : public gaussians {
public:

    gaussians_205(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_205_x : public gaussians {
public:

    dell_gaussians_205_x(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_205_y : public gaussians {
public:

    dell_gaussians_205_y(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class dell_gaussians_205_z : public gaussians {
public:

    dell_gaussians_205_z(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

class lapl_gaussians_205 : public gaussians {
public:

    lapl_gaussians_205(double* exp_factor, double a) : gaussians(exp_factor, a) {
        
    }
        
    virtual double eval(const Walker* walker, int i);

};

/*
    -------------------------  END 99  -------------------------
*/


#endif /* GAUSSIANS_H */
