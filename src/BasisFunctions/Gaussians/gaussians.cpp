
#include "gaussians.h"
#include "../../Walker/Walker.h"


//Superclass Constructor
gaussians::gaussians(double* exp_factor, double a) {
    this->exp_factor = exp_factor;
    this->a = a;
}


/*
    Subclass Eval functions
*/


double gaussians_000::eval(const Walker* walker, int i) {

    (void) walker;
    (void) i;

    
    //exp(-a*r^2)
    
    P = 1;
    return P*(*exp_factor);
    
}

double dell_gaussians_000_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //-2*a*x*exp(-a*r^2)
    
    P = -2*a*x;
    return P*(*exp_factor);
    
}

double dell_gaussians_000_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //-2*a*y*exp(-a*r^2)
    
    P = -2*a*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_000_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //-2*a*z*exp(-a*r^2)
    
    P = -2*a*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_000::eval(const Walker* walker, int i) {
    
    //2*a*(2*a*r^2 - 3)*exp(-a*r^2)
    
    P = 2*a*(2*a*walker->get_r_i2(i) - 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 0  -------------------------
*/

double gaussians_001::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //z*exp(-a*r^2)
    
    P = z;
    return P*(*exp_factor);
    
}

double dell_gaussians_001_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-2*a*x*z*exp(-a*r^2)
    
    P = -2*a*x*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_001_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-2*a*y*z*exp(-a*r^2)
    
    P = -2*a*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_001_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = -2*a*z2 + 1;
    return P*(*exp_factor);
    
}

double lapl_gaussians_001::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //2*a*z*(2*a*r^2 - 5)*exp(-a*r^2)
    
    P = 2*a*z*(2*a*walker->get_r_i2(i) - 5);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 1  -------------------------
*/

double gaussians_010::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //y*exp(-a*r^2)
    
    P = y;
    return P*(*exp_factor);
    
}

double dell_gaussians_010_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-2*a*x*y*exp(-a*r^2)
    
    P = -2*a*x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_010_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = -2*a*y2 + 1;
    return P*(*exp_factor);
    
}

double dell_gaussians_010_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-2*a*y*z*exp(-a*r^2)
    
    P = -2*a*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_010::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //2*a*y*(2*a*r^2 - 5)*exp(-a*r^2)
    
    P = 2*a*y*(2*a*walker->get_r_i2(i) - 5);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 2  -------------------------
*/

double gaussians_100::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //x*exp(-a*r^2)
    
    P = x;
    return P*(*exp_factor);
    
}

double dell_gaussians_100_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = -2*a*x2 + 1;
    return P*(*exp_factor);
    
}

double dell_gaussians_100_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-2*a*x*y*exp(-a*r^2)
    
    P = -2*a*x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_100_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-2*a*x*z*exp(-a*r^2)
    
    P = -2*a*x*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_100::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //2*a*x*(2*a*r^2 - 5)*exp(-a*r^2)
    
    P = 2*a*x*(2*a*walker->get_r_i2(i) - 5);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 3  -------------------------
*/

double gaussians_002::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^2*exp(-a*r^2)
    
    P = z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_002_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*z^2*exp(-a*r^2)
    
    P = -2*a*x*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_002_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*y*z^2*exp(-a*r^2)
    
    P = -2*a*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_002_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_002::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //(z^2*(4*a^2*r^2 - 14*a) + 2)*exp(-a*r^2)
    
    P = z2*(4*pow(a, 2)*walker->get_r_i2(i) - 14*a) + 2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 4  -------------------------
*/

