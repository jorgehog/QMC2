
        
#ifndef HYDROGENIC_H
#define HYDROGENIC_H 

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


#endif /* HYDROGENIC_H */
