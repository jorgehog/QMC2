
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

double gaussians_011::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //y*z*exp(-a*r^2)
    
    P = y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_011_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-2*a*x*y*z*exp(-a*r^2)
    
    P = -2*a*x*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_011_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //z*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_011_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_011::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //2*a*y*z*(2*a*r^2 - 7)*exp(-a*r^2)
    
    P = 2*a*y*z*(2*a*walker->get_r_i2(i) - 7);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 5  -------------------------
*/

double gaussians_020::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^2*exp(-a*r^2)
    
    P = y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_020_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-2*a*x*y^2*exp(-a*r^2)
    
    P = -2*a*x*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_020_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*y*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_020_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*y^2*z*exp(-a*r^2)
    
    P = -2*a*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_020::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //(y^2*(4*a^2*r^2 - 14*a) + 2)*exp(-a*r^2)
    
    P = y2*(4*pow(a, 2)*walker->get_r_i2(i) - 14*a) + 2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 6  -------------------------
*/

double gaussians_101::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //x*z*exp(-a*r^2)
    
    P = x*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_101_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //z*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_101_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-2*a*x*y*z*exp(-a*r^2)
    
    P = -2*a*x*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_101_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_101::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //2*a*x*z*(2*a*r^2 - 7)*exp(-a*r^2)
    
    P = 2*a*x*z*(2*a*walker->get_r_i2(i) - 7);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 7  -------------------------
*/

double gaussians_110::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //x*y*exp(-a*r^2)
    
    P = x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_110_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //y*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_110_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_110_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-2*a*x*y*z*exp(-a*r^2)
    
    P = -2*a*x*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_110::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //2*a*x*y*(2*a*r^2 - 7)*exp(-a*r^2)
    
    P = 2*a*x*y*(2*a*walker->get_r_i2(i) - 7);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 8  -------------------------
*/

double gaussians_200::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^2*exp(-a*r^2)
    
    P = x2;
    return P*(*exp_factor);
    
}

double dell_gaussians_200_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_200_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //-2*a*x^2*y*exp(-a*r^2)
    
    P = -2*a*x2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_200_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^2*z*exp(-a*r^2)
    
    P = -2*a*x2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_200::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //(x^2*(4*a^2*r^2 - 14*a) + 2)*exp(-a*r^2)
    
    P = x2*(4*pow(a, 2)*walker->get_r_i2(i) - 14*a) + 2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 9  -------------------------
*/

double gaussians_003::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^3*exp(-a*r^2)
    
    P = z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_003_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*z^3*exp(-a*r^2)
    
    P = -2*a*x*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_003_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*y*z^3*exp(-a*r^2)
    
    P = -2*a*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_003_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_003::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z*(z^2*(2*a^2*r^2 - 9*a) + 3)*exp(-a*r^2)
    
    P = 2*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 10  -------------------------
*/

double gaussians_012::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^2*exp(-a*r^2)
    
    P = y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_012_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^2*exp(-a*r^2)
    
    P = -2*a*x*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_012_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //z^2*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_012_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_012::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*(z^2*(2*a^2*r^2 - 9*a) + 1)*exp(-a*r^2)
    
    P = 2*y*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 11  -------------------------
*/

double gaussians_021::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^2*z*exp(-a*r^2)
    
    P = y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_021_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^2*z*exp(-a*r^2)
    
    P = -2*a*x*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_021_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y*z*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*y*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_021_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_021::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*z*(y^2*(2*a^2*r^2 - 9*a) + 1)*exp(-a*r^2)
    
    P = 2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 12  -------------------------
*/

double gaussians_030::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^3*exp(-a*r^2)
    
    P = y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_030_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-2*a*x*y^3*exp(-a*r^2)
    
    P = -2*a*x*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_030_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^2*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = y2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_030_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*y^3*z*exp(-a*r^2)
    
    P = -2*a*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_030::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y*(y^2*(2*a^2*r^2 - 9*a) + 3)*exp(-a*r^2)
    
    P = 2*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 13  -------------------------
*/

double gaussians_102::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^2*exp(-a*r^2)
    
    P = x*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_102_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //z^2*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_102_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^2*exp(-a*r^2)
    
    P = -2*a*x*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_102_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_102::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*(z^2*(2*a^2*r^2 - 9*a) + 1)*exp(-a*r^2)
    
    P = 2*x*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 14  -------------------------
*/

double gaussians_111::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //x*y*z*exp(-a*r^2)
    
    P = x*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_111_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //y*z*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_111_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*z*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_111_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_111::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //2*a*x*y*z*(2*a*r^2 - 9)*exp(-a*r^2)
    
    P = 2*a*x*y*z*(2*a*walker->get_r_i2(i) - 9);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 15  -------------------------
*/

double gaussians_120::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^2*exp(-a*r^2)
    
    P = x*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_120_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //y^2*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_120_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_120_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^2*z*exp(-a*r^2)
    
    P = -2*a*x*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_120::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*(y^2*(2*a^2*r^2 - 9*a) + 1)*exp(-a*r^2)
    
    P = 2*x*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 16  -------------------------
*/

double gaussians_201::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^2*z*exp(-a*r^2)
    
    P = x2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_201_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x*z*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_201_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^2*y*z*exp(-a*r^2)
    
    P = -2*a*x2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_201_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_201::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*z*(x^2*(2*a^2*r^2 - 9*a) + 1)*exp(-a*r^2)
    
    P = 2*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 17  -------------------------
*/

double gaussians_210::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //x^2*y*exp(-a*r^2)
    
    P = x2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_210_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //2*x*y*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_210_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_210_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^2*y*z*exp(-a*r^2)
    
    P = -2*a*x2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_210::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //2*y*(x^2*(2*a^2*r^2 - 9*a) + 1)*exp(-a*r^2)
    
    P = 2*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 18  -------------------------
*/

double gaussians_300::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^3*exp(-a*r^2)
    
    P = x2*x;
    return P*(*exp_factor);
    
}

double dell_gaussians_300_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^2*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_300_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //-2*a*x^3*y*exp(-a*r^2)
    
    P = -2*a*x2*x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_300_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^3*z*exp(-a*r^2)
    
    P = -2*a*x2*x*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_300::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x*(x^2*(2*a^2*r^2 - 9*a) + 3)*exp(-a*r^2)
    
    P = 2*x*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 9*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 19  -------------------------
*/

double gaussians_004::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^4*exp(-a*r^2)
    
    P = z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_004_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*z^4*exp(-a*r^2)
    
    P = -2*a*x*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_004_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*y*z^4*exp(-a*r^2)
    
    P = -2*a*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_004_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_004::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z^2*(z^2*(2*a^2*r^2 - 11*a) + 6)*exp(-a*r^2)
    
    P = 2*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 20  -------------------------
*/

double gaussians_013::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^3*exp(-a*r^2)
    
    P = y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_013_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^3*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_013_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //z^3*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = z2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_013_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = y*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_013::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z*(z^2*(2*a^2*r^2 - 11*a) + 3)*exp(-a*r^2)
    
    P = 2*y*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 21  -------------------------
*/

double gaussians_022::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^2*exp(-a*r^2)
    
    P = y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_022_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^2*z^2*exp(-a*r^2)
    
    P = -2*a*x*y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_022_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*z^2*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*y*z2*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_022_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^2*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*y2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_022::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //(y^2*(4*a^2*r^2*z^2 - 22*a*z^2 + 2) + 2*z^2)*exp(-a*r^2)
    
    P = y2*(4*pow(a, 2)*walker->get_r_i2(i)*z2 - 22*a*z2 + 2) + 2*z2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 22  -------------------------
*/

double gaussians_031::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^3*z*exp(-a*r^2)
    
    P = y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_031_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^3*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_031_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^2*z*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = y2*z*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_031_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^3*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = y2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_031::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y*z*(y^2*(2*a^2*r^2 - 11*a) + 3)*exp(-a*r^2)
    
    P = 2*y*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 23  -------------------------
*/

double gaussians_040::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^4*exp(-a*r^2)
    
    P = y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_040_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-2*a*x*y^4*exp(-a*r^2)
    
    P = -2*a*x*y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_040_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y^3*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*y2*y*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_040_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*y^4*z*exp(-a*r^2)
    
    P = -2*a*y2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_040::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y^2*(y^2*(2*a^2*r^2 - 11*a) + 6)*exp(-a*r^2)
    
    P = 2*y2*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 24  -------------------------
*/

double gaussians_103::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^3*exp(-a*r^2)
    
    P = x*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_103_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //z^3*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = z2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_103_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^3*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_103_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_103::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z*(z^2*(2*a^2*r^2 - 11*a) + 3)*exp(-a*r^2)
    
    P = 2*x*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 25  -------------------------
*/

double gaussians_112::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*z^2*exp(-a*r^2)
    
    P = x*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_112_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //y*z^2*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_112_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*z^2*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x*z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_112_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*y*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_112::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*y*(z^2*(2*a^2*r^2 - 11*a) + 1)*exp(-a*r^2)
    
    P = 2*x*y*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 26  -------------------------
*/

double gaussians_121::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*y^2*z*exp(-a*r^2)
    
    P = x*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_121_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //y^2*z*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_121_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*x*y*z*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_121_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x*y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_121::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*x*z*(y^2*(2*a^2*r^2 - 11*a) + 1)*exp(-a*r^2)
    
    P = 2*x*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 27  -------------------------
*/

double gaussians_130::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^3*exp(-a*r^2)
    
    P = x*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_130_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //y^3*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_130_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^2*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x*y2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_130_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^3*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_130::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y*(y^2*(2*a^2*r^2 - 11*a) + 3)*exp(-a*r^2)
    
    P = 2*x*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 28  -------------------------
*/

double gaussians_202::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^2*exp(-a*r^2)
    
    P = x2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_202_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*z^2*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*z2*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_202_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^2*y*z^2*exp(-a*r^2)
    
    P = -2*a*x2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_202_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^2*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_202::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //(x^2*(4*a^2*r^2*z^2 - 22*a*z^2 + 2) + 2*z^2)*exp(-a*r^2)
    
    P = x2*(4*pow(a, 2)*walker->get_r_i2(i)*z2 - 22*a*z2 + 2) + 2*z2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 29  -------------------------
*/

double gaussians_211::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^2*y*z*exp(-a*r^2)
    
    P = x2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_211_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x*y*z*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_211_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*z*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_211_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*y*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_211::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*y*z*(x^2*(2*a^2*r^2 - 11*a) + 1)*exp(-a*r^2)
    
    P = 2*y*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 1);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 30  -------------------------
*/

double gaussians_220::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^2*exp(-a*r^2)
    
    P = x2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_220_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*y^2*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_220_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^2*y*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*y*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_220_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-2*a*x^2*y^2*z*exp(-a*r^2)
    
    P = -2*a*x2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_220::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //(x^2*(4*a^2*r^2*y^2 - 22*a*y^2 + 2) + 2*y^2)*exp(-a*r^2)
    
    P = x2*(4*pow(a, 2)*walker->get_r_i2(i)*y2 - 22*a*y2 + 2) + 2*y2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 31  -------------------------
*/

double gaussians_301::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^3*z*exp(-a*r^2)
    
    P = x2*x*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_301_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^2*z*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*z*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_301_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^3*y*z*exp(-a*r^2)
    
    P = -2*a*x2*x*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_301_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^3*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*x*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_301::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x*z*(x^2*(2*a^2*r^2 - 11*a) + 3)*exp(-a*r^2)
    
    P = 2*x*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 32  -------------------------
*/

double gaussians_310::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //x^3*y*exp(-a*r^2)
    
    P = x2*x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_310_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //x^2*y*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*y*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_310_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^3*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*x*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_310_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^3*y*z*exp(-a*r^2)
    
    P = -2*a*x2*x*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_310::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //2*x*y*(x^2*(2*a^2*r^2 - 11*a) + 3)*exp(-a*r^2)
    
    P = 2*x*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 33  -------------------------
*/

double gaussians_400::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^4*exp(-a*r^2)
    
    P = x2*x2;
    return P*(*exp_factor);
    
}

double dell_gaussians_400_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x^3*(-a*x^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*x*(-a*x2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_400_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //-2*a*x^4*y*exp(-a*r^2)
    
    P = -2*a*x2*x2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_400_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^4*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_400::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x^2*(x^2*(2*a^2*r^2 - 11*a) + 6)*exp(-a*r^2)
    
    P = 2*x2*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 11*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 34  -------------------------
*/

double gaussians_005::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^5*exp(-a*r^2)
    
    P = z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_005_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*z^5*exp(-a*r^2)
    
    P = -2*a*x*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_005_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*y*z^5*exp(-a*r^2)
    
    P = -2*a*y*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_005_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^4*(-2*a*z^2 + 5)*exp(-a*r^2)
    
    P = z2*z2*(-2*a*z2 + 5);
    return P*(*exp_factor);
    
}

double lapl_gaussians_005::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z^3*(z^2*(2*a^2*r^2 - 13*a) + 10)*exp(-a*r^2)
    
    P = 2*z2*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 35  -------------------------
*/

double gaussians_014::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^4*exp(-a*r^2)
    
    P = y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_014_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^4*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_014_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //z^4*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = z2*z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_014_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*y*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_014::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z^2*(z^2*(2*a^2*r^2 - 13*a) + 6)*exp(-a*r^2)
    
    P = 2*y*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 36  -------------------------
*/

double gaussians_023::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^3*exp(-a*r^2)
    
    P = y2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_023_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^2*z^3*exp(-a*r^2)
    
    P = -2*a*x*y2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_023_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*z^3*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*y*z2*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_023_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = y2*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_023::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*z*(y^2*(2*a^2*r^2*z^2 - 13*a*z^2 + 3) + z^2)*exp(-a*r^2)
    
    P = 2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 13*a*z2 + 3) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 37  -------------------------
*/

double gaussians_032::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^3*z^2*exp(-a*r^2)
    
    P = y2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_032_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^3*z^2*exp(-a*r^2)
    
    P = -2*a*x*y2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_032_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^2*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = y2*z2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_032_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^3*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*y2*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_032::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*(y^2*(2*a^2*r^2*z^2 - 13*a*z^2 + 1) + 3*z^2)*exp(-a*r^2)
    
    P = 2*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 13*a*z2 + 1) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 38  -------------------------
*/

double gaussians_041::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^4*z*exp(-a*r^2)
    
    P = y2*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_041_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^4*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_041_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y^3*z*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*y2*y*z*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_041_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^4*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_041::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y^2*z*(y^2*(2*a^2*r^2 - 13*a) + 6)*exp(-a*r^2)
    
    P = 2*y2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 39  -------------------------
*/

double gaussians_050::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^5*exp(-a*r^2)
    
    P = y2*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_050_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-2*a*x*y^5*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_050_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^4*(-2*a*y^2 + 5)*exp(-a*r^2)
    
    P = y2*y2*(-2*a*y2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_050_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*y^5*z*exp(-a*r^2)
    
    P = -2*a*y2*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_050::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y^3*(y^2*(2*a^2*r^2 - 13*a) + 10)*exp(-a*r^2)
    
    P = 2*y2*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 40  -------------------------
*/

double gaussians_104::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^4*exp(-a*r^2)
    
    P = x*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_104_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //z^4*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = z2*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_104_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^4*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_104_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*x*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_104::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z^2*(z^2*(2*a^2*r^2 - 13*a) + 6)*exp(-a*r^2)
    
    P = 2*x*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 41  -------------------------
*/

double gaussians_113::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*z^3*exp(-a*r^2)
    
    P = x*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_113_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //y*z^3*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y*z2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_113_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*z^3*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x*z2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_113_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x*y*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_113::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*y*z*(z^2*(2*a^2*r^2 - 13*a) + 3)*exp(-a*r^2)
    
    P = 2*x*y*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 42  -------------------------
*/

double gaussians_122::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*z^2*exp(-a*r^2)
    
    P = x*y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_122_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^2*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_122_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y*z^2*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z2*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_122_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^2*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_122::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*(y^2*(2*a^2*r^2*z^2 - 13*a*z^2 + 1) + z^2)*exp(-a*r^2)
    
    P = 2*x*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 13*a*z2 + 1) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 43  -------------------------
*/

double gaussians_131::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*y^3*z*exp(-a*r^2)
    
    P = x*y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_131_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //y^3*z*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_131_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*y^2*z*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x*y2*z*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_131_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^3*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x*y2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_131::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*x*y*z*(y^2*(2*a^2*r^2 - 13*a) + 3)*exp(-a*r^2)
    
    P = 2*x*y*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 44  -------------------------
*/

double gaussians_140::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^4*exp(-a*r^2)
    
    P = x*y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_140_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //y^4*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_140_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y^3*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*x*y2*y*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_140_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^4*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_140::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y^2*(y^2*(2*a^2*r^2 - 13*a) + 6)*exp(-a*r^2)
    
    P = 2*x*y2*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 45  -------------------------
*/

double gaussians_203::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^3*exp(-a*r^2)
    
    P = x2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_203_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*z^3*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*z2*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_203_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^2*y*z^3*exp(-a*r^2)
    
    P = -2*a*x2*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_203_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x2*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_203::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*z*(x^2*(2*a^2*r^2*z^2 - 13*a*z^2 + 3) + z^2)*exp(-a*r^2)
    
    P = 2*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 13*a*z2 + 3) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 46  -------------------------
*/

double gaussians_212::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*y*z^2*exp(-a*r^2)
    
    P = x2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_212_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*y*z^2*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z2*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_212_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^2*z^2*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_212_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^2*y*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_212::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*y*(x^2*(2*a^2*r^2*z^2 - 13*a*z^2 + 1) + z^2)*exp(-a*r^2)
    
    P = 2*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 13*a*z2 + 1) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 47  -------------------------
*/

double gaussians_221::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^2*z*exp(-a*r^2)
    
    P = x2*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_221_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*y^2*z*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_221_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^2*y*z*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*y*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_221_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^2*y^2*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_221::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*z*(x^2*(2*a^2*r^2*y^2 - 13*a*y^2 + 1) + y^2)*exp(-a*r^2)
    
    P = 2*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 13*a*y2 + 1) + y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 48  -------------------------
*/

double gaussians_230::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^3*exp(-a*r^2)
    
    P = x2*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_230_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*y^3*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*y*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_230_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^2*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x2*y2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_230_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-2*a*x^2*y^3*z*exp(-a*r^2)
    
    P = -2*a*x2*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_230::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*y*(x^2*(2*a^2*r^2*y^2 - 13*a*y^2 + 3) + y^2)*exp(-a*r^2)
    
    P = 2*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 13*a*y2 + 3) + y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 49  -------------------------
*/

double gaussians_302::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^3*z^2*exp(-a*r^2)
    
    P = x2*x*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_302_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^2*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*z2*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_302_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^3*y*z^2*exp(-a*r^2)
    
    P = -2*a*x2*x*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_302_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^3*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*x*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_302::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*(x^2*(2*a^2*r^2*z^2 - 13*a*z^2 + 1) + 3*z^2)*exp(-a*r^2)
    
    P = 2*x*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 13*a*z2 + 1) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 50  -------------------------
*/

double gaussians_311::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^3*y*z*exp(-a*r^2)
    
    P = x2*x*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_311_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^2*y*z*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*y*z*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_311_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^3*z*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*x*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_311_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^3*y*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*x*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_311::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x*y*z*(x^2*(2*a^2*r^2 - 13*a) + 3)*exp(-a*r^2)
    
    P = 2*x*y*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 3);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 51  -------------------------
*/

double gaussians_320::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^3*y^2*exp(-a*r^2)
    
    P = x2*x*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_320_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^2*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*y2*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_320_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^3*y*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*x*y*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_320_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-2*a*x^3*y^2*z*exp(-a*r^2)
    
    P = -2*a*x2*x*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_320::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*(x^2*(2*a^2*r^2*y^2 - 13*a*y^2 + 1) + 3*y^2)*exp(-a*r^2)
    
    P = 2*x*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 13*a*y2 + 1) + 3*y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 52  -------------------------
*/

double gaussians_401::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^4*z*exp(-a*r^2)
    
    P = x2*x2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_401_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x^3*z*(-a*x^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*x*z*(-a*x2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_401_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^4*y*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_401_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^4*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*x2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_401::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x^2*z*(x^2*(2*a^2*r^2 - 13*a) + 6)*exp(-a*r^2)
    
    P = 2*x2*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 53  -------------------------
*/

double gaussians_410::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //x^4*y*exp(-a*r^2)
    
    P = x2*x2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_410_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //2*x^3*y*(-a*x^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*x*y*(-a*x2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_410_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^4*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*x2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_410_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^4*y*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_410::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //2*x^2*y*(x^2*(2*a^2*r^2 - 13*a) + 6)*exp(-a*r^2)
    
    P = 2*x2*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 54  -------------------------
*/

double gaussians_500::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^5*exp(-a*r^2)
    
    P = x2*x2*x;
    return P*(*exp_factor);
    
}

double dell_gaussians_500_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^4*(-2*a*x^2 + 5)*exp(-a*r^2)
    
    P = x2*x2*(-2*a*x2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_500_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //-2*a*x^5*y*exp(-a*r^2)
    
    P = -2*a*x2*x2*x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_500_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^5*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*x*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_500::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x^3*(x^2*(2*a^2*r^2 - 13*a) + 10)*exp(-a*r^2)
    
    P = 2*x2*x*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 13*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 55  -------------------------
*/

double gaussians_006::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^6*exp(-a*r^2)
    
    P = z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_006_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*z^6*exp(-a*r^2)
    
    P = -2*a*x*z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_006_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*y*z^6*exp(-a*r^2)
    
    P = -2*a*y*z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_006_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z^5*(-a*z^2 + 3)*exp(-a*r^2)
    
    P = 2*z2*z2*z*(-a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_006::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z^4*(z^2*(2*a^2*r^2 - 15*a) + 15)*exp(-a*r^2)
    
    P = 2*z2*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 56  -------------------------
*/

double gaussians_015::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^5*exp(-a*r^2)
    
    P = y*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_015_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^5*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_015_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //z^5*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = z2*z2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_015_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^4*(-2*a*z^2 + 5)*exp(-a*r^2)
    
    P = y*z2*z2*(-2*a*z2 + 5);
    return P*(*exp_factor);
    
}

double lapl_gaussians_015::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z^3*(z^2*(2*a^2*r^2 - 15*a) + 10)*exp(-a*r^2)
    
    P = 2*y*z2*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 57  -------------------------
*/

double gaussians_024::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^4*exp(-a*r^2)
    
    P = y2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_024_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^2*z^4*exp(-a*r^2)
    
    P = -2*a*x*y2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_024_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*z^4*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*y*z2*z2*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_024_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^2*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*y2*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_024::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*z^2*(y^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 6) + z^2)*exp(-a*r^2)
    
    P = 2*z2*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 6) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 58  -------------------------
*/

double gaussians_033::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^3*z^3*exp(-a*r^2)
    
    P = y2*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_033_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^3*z^3*exp(-a*r^2)
    
    P = -2*a*x*y2*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_033_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^3*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = y2*z2*z*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_033_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^3*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = y2*y*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_033::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*z*(y^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 3) + 3*z^2)*exp(-a*r^2)
    
    P = 2*y*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 3) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 59  -------------------------
*/

double gaussians_042::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^4*z^2*exp(-a*r^2)
    
    P = y2*y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_042_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^4*z^2*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_042_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^3*z^2*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*y2*y*z2*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_042_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^4*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*y2*y2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_042::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^2*(y^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 1) + 6*z^2)*exp(-a*r^2)
    
    P = 2*y2*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 1) + 6*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 60  -------------------------
*/

double gaussians_051::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^5*z*exp(-a*r^2)
    
    P = y2*y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_051_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^5*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_051_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^4*z*(-2*a*y^2 + 5)*exp(-a*r^2)
    
    P = y2*y2*z*(-2*a*y2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_051_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^5*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_051::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y^3*z*(y^2*(2*a^2*r^2 - 15*a) + 10)*exp(-a*r^2)
    
    P = 2*y2*y*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 61  -------------------------
*/

double gaussians_060::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^6*exp(-a*r^2)
    
    P = y2*y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_060_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-2*a*x*y^6*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_060_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y^5*(-a*y^2 + 3)*exp(-a*r^2)
    
    P = 2*y2*y2*y*(-a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_060_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*y^6*z*exp(-a*r^2)
    
    P = -2*a*y2*y2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_060::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y^4*(y^2*(2*a^2*r^2 - 15*a) + 15)*exp(-a*r^2)
    
    P = 2*y2*y2*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 62  -------------------------
*/

double gaussians_105::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^5*exp(-a*r^2)
    
    P = x*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_105_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //z^5*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = z2*z2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_105_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^5*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_105_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^4*(-2*a*z^2 + 5)*exp(-a*r^2)
    
    P = x*z2*z2*(-2*a*z2 + 5);
    return P*(*exp_factor);
    
}

double lapl_gaussians_105::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z^3*(z^2*(2*a^2*r^2 - 15*a) + 10)*exp(-a*r^2)
    
    P = 2*x*z2*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 63  -------------------------
*/

double gaussians_114::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*z^4*exp(-a*r^2)
    
    P = x*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_114_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //y*z^4*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y*z2*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_114_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*z^4*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x*z2*z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_114_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*y*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*x*y*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_114::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*y*z^2*(z^2*(2*a^2*r^2 - 15*a) + 6)*exp(-a*r^2)
    
    P = 2*x*y*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 64  -------------------------
*/

double gaussians_123::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*z^3*exp(-a*r^2)
    
    P = x*y2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_123_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^3*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*z2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_123_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y*z^3*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z2*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_123_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x*y2*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_123::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*z*(y^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 3) + z^2)*exp(-a*r^2)
    
    P = 2*x*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 3) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 65  -------------------------
*/

double gaussians_132::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^3*z^2*exp(-a*r^2)
    
    P = x*y2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_132_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //y^3*z^2*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_132_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*z^2*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x*y2*z2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_132_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^3*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_132::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y*(y^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 1) + 3*z^2)*exp(-a*r^2)
    
    P = 2*x*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 1) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 66  -------------------------
*/

double gaussians_141::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*y^4*z*exp(-a*r^2)
    
    P = x*y2*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_141_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //y^4*z*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_141_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*x*y^3*z*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*x*y2*y*z*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_141_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^4*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x*y2*y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_141::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*x*y^2*z*(y^2*(2*a^2*r^2 - 15*a) + 6)*exp(-a*r^2)
    
    P = 2*x*y2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 67  -------------------------
*/

double gaussians_150::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^5*exp(-a*r^2)
    
    P = x*y2*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_150_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //y^5*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*y*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_150_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^4*(-2*a*y^2 + 5)*exp(-a*r^2)
    
    P = x*y2*y2*(-2*a*y2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_150_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^5*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_150::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y^3*(y^2*(2*a^2*r^2 - 15*a) + 10)*exp(-a*r^2)
    
    P = 2*x*y2*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 68  -------------------------
*/

double gaussians_204::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^4*exp(-a*r^2)
    
    P = x2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_204_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*z^4*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*z2*z2*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_204_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^2*y*z^4*exp(-a*r^2)
    
    P = -2*a*x2*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_204_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^2*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_204::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*z^2*(x^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 6) + z^2)*exp(-a*r^2)
    
    P = 2*z2*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 6) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 69  -------------------------
*/

double gaussians_213::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*y*z^3*exp(-a*r^2)
    
    P = x2*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_213_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*y*z^3*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z2*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_213_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^2*z^3*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*z2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_213_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*y*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x2*y*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_213::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*y*z*(x^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 3) + z^2)*exp(-a*r^2)
    
    P = 2*y*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 3) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 70  -------------------------
*/

double gaussians_222::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^2*y^2*z^2*exp(-a*r^2)
    
    P = x2*y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_222_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^2*z^2*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*z2*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_222_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //2*x^2*y*z^2*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*y*z2*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_222_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //2*x^2*y^2*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*y2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_222::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //(x^2*(4*a^2*r^2*y^2*z^2 - 30*a*y^2*z^2 + 2*y^2 + 2*z^2) + 2*y^2*z^2)*exp(-a*r^2)
    
    P = x2*(4*pow(a, 2)*walker->get_r_i2(i)*y2*z2 - 30*a*y2*z2 + 2*y2 + 2*z2) + 2*y2*z2;
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 71  -------------------------
*/

double gaussians_231::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^3*z*exp(-a*r^2)
    
    P = x2*y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_231_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*y^3*z*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*y*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_231_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^2*z*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x2*y2*z*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_231_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^2*y^3*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*y2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_231::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*y*z*(x^2*(2*a^2*r^2*y^2 - 15*a*y^2 + 3) + y^2)*exp(-a*r^2)
    
    P = 2*y*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 15*a*y2 + 3) + y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 72  -------------------------
*/

double gaussians_240::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^4*exp(-a*r^2)
    
    P = x2*y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_240_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*y^4*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*y2*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_240_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^2*y^3*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*y2*y*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_240_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-2*a*x^2*y^4*z*exp(-a*r^2)
    
    P = -2*a*x2*y2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_240::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*y^2*(x^2*(2*a^2*r^2*y^2 - 15*a*y^2 + 6) + y^2)*exp(-a*r^2)
    
    P = 2*y2*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 15*a*y2 + 6) + y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 73  -------------------------
*/

double gaussians_303::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^3*z^3*exp(-a*r^2)
    
    P = x2*x*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_303_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^3*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*z2*z*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_303_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^3*y*z^3*exp(-a*r^2)
    
    P = -2*a*x2*x*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_303_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^3*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x2*x*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_303::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*z*(x^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 3) + 3*z^2)*exp(-a*r^2)
    
    P = 2*x*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 3) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 74  -------------------------
*/

double gaussians_312::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^3*y*z^2*exp(-a*r^2)
    
    P = x2*x*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_312_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*y*z^2*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*y*z2*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_312_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^3*z^2*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*x*z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_312_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^3*y*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*x*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_312::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*y*(x^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 1) + 3*z^2)*exp(-a*r^2)
    
    P = 2*x*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 1) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 75  -------------------------
*/

double gaussians_321::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^3*y^2*z*exp(-a*r^2)
    
    P = x2*x*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_321_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^2*z*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*y2*z*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_321_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^3*y*z*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*x*y*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_321_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //x^3*y^2*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*x*y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_321::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*z*(x^2*(2*a^2*r^2*y^2 - 15*a*y^2 + 1) + 3*y^2)*exp(-a*r^2)
    
    P = 2*x*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 15*a*y2 + 1) + 3*y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 76  -------------------------
*/

double gaussians_330::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^3*y^3*exp(-a*r^2)
    
    P = x2*x*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_330_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^2*y^3*(-2*a*x^2 + 3)*exp(-a*r^2)
    
    P = x2*y2*y*(-2*a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_330_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^3*y^2*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x2*x*y2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_330_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-2*a*x^3*y^3*z*exp(-a*r^2)
    
    P = -2*a*x2*x*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_330::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x*y*(x^2*(2*a^2*r^2*y^2 - 15*a*y^2 + 3) + 3*y^2)*exp(-a*r^2)
    
    P = 2*x*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 15*a*y2 + 3) + 3*y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 77  -------------------------
*/

double gaussians_402::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^4*z^2*exp(-a*r^2)
    
    P = x2*x2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_402_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^3*z^2*(-a*x^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*x*z2*(-a*x2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_402_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^4*y*z^2*exp(-a*r^2)
    
    P = -2*a*x2*x2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_402_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^4*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*x2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_402::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x^2*(x^2*(2*a^2*r^2*z^2 - 15*a*z^2 + 1) + 6*z^2)*exp(-a*r^2)
    
    P = 2*x2*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 15*a*z2 + 1) + 6*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 78  -------------------------
*/

double gaussians_411::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^4*y*z*exp(-a*r^2)
    
    P = x2*x2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_411_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x^3*y*z*(-a*x^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*x*y*z*(-a*x2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_411_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //x^4*z*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*x2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_411_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^4*y*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*x2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_411::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x^2*y*z*(x^2*(2*a^2*r^2 - 15*a) + 6)*exp(-a*r^2)
    
    P = 2*x2*y*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 6);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 79  -------------------------
*/

double gaussians_420::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^4*y^2*exp(-a*r^2)
    
    P = x2*x2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_420_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^3*y^2*(-a*x^2 + 2)*exp(-a*r^2)
    
    P = 2*x2*x*y2*(-a*x2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_420_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^4*y*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x2*x2*y*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_420_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-2*a*x^4*y^2*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_420::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //2*x^2*(x^2*(2*a^2*r^2*y^2 - 15*a*y^2 + 1) + 6*y^2)*exp(-a*r^2)
    
    P = 2*x2*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*y2 - 15*a*y2 + 1) + 6*y2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 80  -------------------------
*/

double gaussians_501::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^5*z*exp(-a*r^2)
    
    P = x2*x2*x*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_501_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //x^4*z*(-2*a*x^2 + 5)*exp(-a*r^2)
    
    P = x2*x2*z*(-2*a*x2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_501_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^5*y*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*x*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_501_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^5*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x2*x2*x*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_501::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //2*x^3*z*(x^2*(2*a^2*r^2 - 15*a) + 10)*exp(-a*r^2)
    
    P = 2*x2*x*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 81  -------------------------
*/

double gaussians_510::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //x^5*y*exp(-a*r^2)
    
    P = x2*x2*x*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_510_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //x^4*y*(-2*a*x^2 + 5)*exp(-a*r^2)
    
    P = x2*x2*y*(-2*a*x2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_510_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //x^5*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x2*x2*x*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_510_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^5*y*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*x*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_510::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //2*x^3*y*(x^2*(2*a^2*r^2 - 15*a) + 10)*exp(-a*r^2)
    
    P = 2*x2*x*y*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 82  -------------------------
*/

double gaussians_600::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //x^6*exp(-a*r^2)
    
    P = x2*x2*x2;
    return P*(*exp_factor);
    
}

double dell_gaussians_600_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x^5*(-a*x^2 + 3)*exp(-a*r^2)
    
    P = 2*x2*x2*x*(-a*x2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_600_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //-2*a*x^6*y*exp(-a*r^2)
    
    P = -2*a*x2*x2*x2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_600_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-2*a*x^6*z*exp(-a*r^2)
    
    P = -2*a*x2*x2*x2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_600::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //2*x^4*(x^2*(2*a^2*r^2 - 15*a) + 15)*exp(-a*r^2)
    
    P = 2*x2*x2*(x2*(2*pow(a, 2)*walker->get_r_i2(i) - 15*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 83  -------------------------
*/

double gaussians_007::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^7*exp(-a*r^2)
    
    P = z2*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_007_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*z^7*exp(-a*r^2)
    
    P = -2*a*x*z2*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_007_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*y*z^7*exp(-a*r^2)
    
    P = -2*a*y*z2*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_007_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //z^6*(-2*a*z^2 + 7)*exp(-a*r^2)
    
    P = z2*z2*z2*(-2*a*z2 + 7);
    return P*(*exp_factor);
    
}

double lapl_gaussians_007::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*z^5*(z^2*(2*a^2*r^2 - 17*a) + 21)*exp(-a*r^2)
    
    P = 2*z2*z2*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 21);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 84  -------------------------
*/

double gaussians_016::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*z^6*exp(-a*r^2)
    
    P = y*z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_016_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^6*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_016_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //z^6*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = z2*z2*z2*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_016_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z^5*(-a*z^2 + 3)*exp(-a*r^2)
    
    P = 2*y*z2*z2*z*(-a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_016::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*y*z^4*(z^2*(2*a^2*r^2 - 17*a) + 15)*exp(-a*r^2)
    
    P = 2*y*z2*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 85  -------------------------
*/

double gaussians_025::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^5*exp(-a*r^2)
    
    P = y2*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_025_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^2*z^5*exp(-a*r^2)
    
    P = -2*a*x*y2*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_025_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*z^5*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*y*z2*z2*z*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_025_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^4*(-2*a*z^2 + 5)*exp(-a*r^2)
    
    P = y2*z2*z2*(-2*a*z2 + 5);
    return P*(*exp_factor);
    
}

double lapl_gaussians_025::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*z^3*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 10) + z^2)*exp(-a*r^2)
    
    P = 2*z2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 10) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 86  -------------------------
*/

double gaussians_034::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^3*z^4*exp(-a*r^2)
    
    P = y2*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_034_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^3*z^4*exp(-a*r^2)
    
    P = -2*a*x*y2*y*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_034_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^4*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = y2*z2*z2*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_034_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^3*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*y2*y*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_034::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y*z^2*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 6) + 3*z^2)*exp(-a*r^2)
    
    P = 2*y*z2*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 6) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 87  -------------------------
*/

double gaussians_043::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^4*z^3*exp(-a*r^2)
    
    P = y2*y2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_043_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^4*z^3*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_043_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^3*z^3*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*y2*y*z2*z*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_043_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^4*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = y2*y2*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_043::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^2*z*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 3) + 6*z^2)*exp(-a*r^2)
    
    P = 2*y2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 3) + 6*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 88  -------------------------
*/

double gaussians_052::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^5*z^2*exp(-a*r^2)
    
    P = y2*y2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_052_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //-2*a*x*y^5*z^2*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_052_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^4*z^2*(-2*a*y^2 + 5)*exp(-a*r^2)
    
    P = y2*y2*z2*(-2*a*y2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_052_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^5*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*y2*y2*y*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_052::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*y^3*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 1) + 10*z^2)*exp(-a*r^2)
    
    P = 2*y2*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 1) + 10*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 89  -------------------------
*/

double gaussians_061::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //y^6*z*exp(-a*r^2)
    
    P = y2*y2*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_061_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^6*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_061_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y^5*z*(-a*y^2 + 3)*exp(-a*r^2)
    
    P = 2*y2*y2*y*z*(-a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_061_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //y^6*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*y2*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_061::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*y^4*z*(y^2*(2*a^2*r^2 - 17*a) + 15)*exp(-a*r^2)
    
    P = 2*y2*y2*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 90  -------------------------
*/

double gaussians_070::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^7*exp(-a*r^2)
    
    P = y2*y2*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_070_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-2*a*x*y^7*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y2*y;
    return P*(*exp_factor);
    
}

double dell_gaussians_070_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //y^6*(-2*a*y^2 + 7)*exp(-a*r^2)
    
    P = y2*y2*y2*(-2*a*y2 + 7);
    return P*(*exp_factor);
    
}

double dell_gaussians_070_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*y^7*z*exp(-a*r^2)
    
    P = -2*a*y2*y2*y2*y*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_070::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*y^5*(y^2*(2*a^2*r^2 - 17*a) + 21)*exp(-a*r^2)
    
    P = 2*y2*y2*y*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 21);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 91  -------------------------
*/

double gaussians_106::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*z^6*exp(-a*r^2)
    
    P = x*z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_106_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //z^6*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = z2*z2*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_106_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-2*a*x*y*z^6*exp(-a*r^2)
    
    P = -2*a*x*y*z2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_106_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z^5*(-a*z^2 + 3)*exp(-a*r^2)
    
    P = 2*x*z2*z2*z*(-a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_106::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*z^4*(z^2*(2*a^2*r^2 - 17*a) + 15)*exp(-a*r^2)
    
    P = 2*x*z2*z2*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 92  -------------------------
*/

double gaussians_115::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*z^5*exp(-a*r^2)
    
    P = x*y*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_115_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //y*z^5*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y*z2*z2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_115_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*z^5*(-2*a*y^2 + 1)*exp(-a*r^2)
    
    P = x*z2*z2*z*(-2*a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_115_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //x*y*z^4*(-2*a*z^2 + 5)*exp(-a*r^2)
    
    P = x*y*z2*z2*(-2*a*z2 + 5);
    return P*(*exp_factor);
    
}

double lapl_gaussians_115::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //2*x*y*z^3*(z^2*(2*a^2*r^2 - 17*a) + 10)*exp(-a*r^2)
    
    P = 2*x*y*z2*z*(z2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 93  -------------------------
*/

double gaussians_124::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*z^4*exp(-a*r^2)
    
    P = x*y2*z2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_124_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //y^2*z^4*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*z2*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_124_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y*z^4*(-a*y^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y*z2*z2*(-a*y2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_124_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^2*z^3*(-a*z^2 + 2)*exp(-a*r^2)
    
    P = 2*x*y2*z2*z*(-a*z2 + 2);
    return P*(*exp_factor);
    
}

double lapl_gaussians_124::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*z^2*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 6) + z^2)*exp(-a*r^2)
    
    P = 2*x*z2*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 6) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 94  -------------------------
*/

double gaussians_133::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^3*z^3*exp(-a*r^2)
    
    P = x*y2*y*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_133_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //y^3*z^3*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y*z2*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_133_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^2*z^3*(-2*a*y^2 + 3)*exp(-a*r^2)
    
    P = x*y2*z2*z*(-2*a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_133_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^3*z^2*(-2*a*z^2 + 3)*exp(-a*r^2)
    
    P = x*y2*y*z2*(-2*a*z2 + 3);
    return P*(*exp_factor);
    
}

double lapl_gaussians_133::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y*z*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 3) + 3*z^2)*exp(-a*r^2)
    
    P = 2*x*y*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 3) + 3*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 95  -------------------------
*/

double gaussians_142::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^4*z^2*exp(-a*r^2)
    
    P = x*y2*y2*z2;
    return P*(*exp_factor);
    
}

double dell_gaussians_142_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    
    //y^4*z^2*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*z2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_142_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^3*z^2*(-a*y^2 + 2)*exp(-a*r^2)
    
    P = 2*x*y2*y*z2*(-a*y2 + 2);
    return P*(*exp_factor);
    
}

double dell_gaussians_142_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^4*z*(-a*z^2 + 1)*exp(-a*r^2)
    
    P = 2*x*y2*y2*z*(-a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_142::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //2*x*y^2*(y^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 1) + 6*z^2)*exp(-a*r^2)
    
    P = 2*x*y2*(y2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 1) + 6*z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 96  -------------------------
*/

double gaussians_151::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*y^5*z*exp(-a*r^2)
    
    P = x*y2*y2*y*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_151_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //y^5*z*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*y*z*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_151_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //x*y^4*z*(-2*a*y^2 + 5)*exp(-a*r^2)
    
    P = x*y2*y2*z*(-2*a*y2 + 5);
    return P*(*exp_factor);
    
}

double dell_gaussians_151_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    z2 = z*z;
    
    //x*y^5*(-2*a*z^2 + 1)*exp(-a*r^2)
    
    P = x*y2*y2*y*(-2*a*z2 + 1);
    return P*(*exp_factor);
    
}

double lapl_gaussians_151::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //2*x*y^3*z*(y^2*(2*a^2*r^2 - 17*a) + 10)*exp(-a*r^2)
    
    P = 2*x*y2*y*z*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 10);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 97  -------------------------
*/

double gaussians_160::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //x*y^6*exp(-a*r^2)
    
    P = x*y2*y2*y2;
    return P*(*exp_factor);
    
}

double dell_gaussians_160_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //y^6*(-2*a*x^2 + 1)*exp(-a*r^2)
    
    P = y2*y2*y2*(-2*a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_160_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y^5*(-a*y^2 + 3)*exp(-a*r^2)
    
    P = 2*x*y2*y2*y*(-a*y2 + 3);
    return P*(*exp_factor);
    
}

double dell_gaussians_160_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-2*a*x*y^6*z*exp(-a*r^2)
    
    P = -2*a*x*y2*y2*y2*z;
    return P*(*exp_factor);
    
}

double lapl_gaussians_160::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //2*x*y^4*(y^2*(2*a^2*r^2 - 17*a) + 15)*exp(-a*r^2)
    
    P = 2*x*y2*y2*(y2*(2*pow(a, 2)*walker->get_r_i2(i) - 17*a) + 15);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 98  -------------------------
*/

double gaussians_205::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^5*exp(-a*r^2)
    
    P = x2*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_205_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*x*z^5*(-a*x^2 + 1)*exp(-a*r^2)
    
    P = 2*x*z2*z2*z*(-a*x2 + 1);
    return P*(*exp_factor);
    
}

double dell_gaussians_205_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //-2*a*x^2*y*z^5*exp(-a*r^2)
    
    P = -2*a*x2*y*z2*z2*z;
    return P*(*exp_factor);
    
}

double dell_gaussians_205_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //x^2*z^4*(-2*a*z^2 + 5)*exp(-a*r^2)
    
    P = x2*z2*z2*(-2*a*z2 + 5);
    return P*(*exp_factor);
    
}

double lapl_gaussians_205::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //2*z^3*(x^2*(2*a^2*r^2*z^2 - 17*a*z^2 + 10) + z^2)*exp(-a*r^2)
    
    P = 2*z2*z*(x2*(2*pow(a, 2)*walker->get_r_i2(i)*z2 - 17*a*z2 + 10) + z2);
    return P*(*exp_factor);
    
}

/*
    -------------------------  END 99  -------------------------
*/

