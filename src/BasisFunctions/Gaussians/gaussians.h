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


#endif /* GAUSSIANS_H */
