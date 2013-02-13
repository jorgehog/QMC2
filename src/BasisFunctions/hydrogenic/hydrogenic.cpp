
#include "../../QMCheaders.h"


//Superclass Constructor
hydrogenic::hydrogenic(double* k, double* k2, double* r22d, double* r2d, double* exp_factor) {
    this->k = k;
    this->k2 = k2;
    this->r22d = r22d;
    this->r2d = r2d;
    this->exp_factor = exp_factor;
}



/*
    Subclass Constructors
*/


hydrogenic_0::hydrogenic_0(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_0_x::dell_hydrogenic_0_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_0_y::dell_hydrogenic_0_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_0_z::dell_hydrogenic_0_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_0::lapl_hydrogenic_0(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 0  -------------------------
*/

hydrogenic_1::hydrogenic_1(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_1_x::dell_hydrogenic_1_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_1_y::dell_hydrogenic_1_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_1_z::dell_hydrogenic_1_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_1::lapl_hydrogenic_1(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 1  -------------------------
*/

hydrogenic_2::hydrogenic_2(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_2_x::dell_hydrogenic_2_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_2_y::dell_hydrogenic_2_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_2_z::dell_hydrogenic_2_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_2::lapl_hydrogenic_2(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 2  -------------------------
*/

hydrogenic_3::hydrogenic_3(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_3_x::dell_hydrogenic_3_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_3_y::dell_hydrogenic_3_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_3_z::dell_hydrogenic_3_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_3::lapl_hydrogenic_3(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 3  -------------------------
*/

hydrogenic_4::hydrogenic_4(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_4_x::dell_hydrogenic_4_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_4_y::dell_hydrogenic_4_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_4_z::dell_hydrogenic_4_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_4::lapl_hydrogenic_4(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 4  -------------------------
*/



/*
    Subclass Eval functions
*/


double hydrogenic_0::eval(const Walker* walker, int i) {


    
    //exp(-k*r)
    
    psi = 1;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_0_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //-k*x*exp(-k*r)/r
    
    psi = -(*k)*x/walker->get_r_i(i);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_0_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //-k*y*exp(-k*r)/r
    
    psi = -(*k)*y/walker->get_r_i(i);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_0_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //-k*z*exp(-k*r)/r
    
    psi = -(*k)*z/walker->get_r_i(i);
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_0::eval(const Walker* walker, int i) {


    
    //k*(k*r - 2)*exp(-k*r)/r
    
    psi = (*k)*((*k)*walker->get_r_i(i) - 2)/walker->get_r_i(i);
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 0  -------------------------
*/

double hydrogenic_1::eval(const Walker* walker, int i) {


    
    //(k*r - 2)*exp(-k*r/2)
    
    psi = (*k)*walker->get_r_i(i) - 2;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_1_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //-k*x*(k*r - 4)*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*x*((*k)*walker->get_r_i(i) - 4)/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_1_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //-k*y*(k*r - 4)*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*y*((*k)*walker->get_r_i(i) - 4)/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_1_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //-k*z*(k*r - 4)*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*z*((*k)*walker->get_r_i(i) - 4)/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_1::eval(const Walker* walker, int i) {


    
    //k*(k*r - 8)*(k*r - 2)*exp(-k*r/2)/(4*r)
    
    psi = (*k)*((*k)*walker->get_r_i(i) - 8)*((*k)*walker->get_r_i(i) - 2)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 1  -------------------------
*/

double hydrogenic_2::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //r*x*exp(-k*r/2)/r_2d
    
    psi = walker->get_r_i(i)*x/(*r2d);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_2_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    z2 = z*z;
    
    //(-r_2d^2*(k*r*x^2 - 2*r^2) - 2*x^2*z^2)*exp(-k*r/2)/(2*r*r_2d^3)
    
    psi = (-(*r22d)*((*k)*walker->get_r_i(i)*x2 - 2*walker->get_r_i2(i)) - 2*x2*z2)/(2*walker->get_r_i(i)*(*r22d)*(*r2d));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_2_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-x*y*(k*r*r_2d^2 + 2*z^2)*exp(-k*r/2)/(2*r*r_2d^3)
    
    psi = -x*y*((*k)*walker->get_r_i(i)*(*r22d) + 2*z2)/(2*walker->get_r_i(i)*(*r22d)*(*r2d));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_2_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-x*z*(k*r - 2)*exp(-k*r/2)/(2*r*r_2d)
    
    psi = -x*z*((*k)*walker->get_r_i(i) - 2)/(2*walker->get_r_i(i)*(*r2d));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_2::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //x*(-4*r^2 + r_2d^2*(k^2*r^2 - 8*k*r + 8))*exp(-k*r/2)/(4*r*r_2d^3)
    
    psi = x*(-4*walker->get_r_i2(i) + (*r22d)*((*k2)*walker->get_r_i2(i) - 8*(*k)*walker->get_r_i(i) + 8))/(4*walker->get_r_i(i)*(*r22d)*(*r2d));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 2  -------------------------
*/

double hydrogenic_3::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //y*z*exp(-k*r/2)/r_2d
    
    psi = y*z/(*r2d);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_3_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-x*y*z*(k*r_2d^2 + 2*r)*exp(-k*r/2)/(2*r*r_2d^3)
    
    psi = -x*y*z*((*k)*(*r22d) + 2*walker->get_r_i(i))/(2*walker->get_r_i(i)*(*r22d)*(*r2d));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_3_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-z*(2*r*y^2 + r_2d^2*(k*y^2 - 2*r))*exp(-k*r/2)/(2*r*r_2d^3)
    
    psi = -z*(2*walker->get_r_i(i)*y2 + (*r22d)*((*k)*y2 - 2*walker->get_r_i(i)))/(2*walker->get_r_i(i)*(*r22d)*(*r2d));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_3_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //y*(-k*z^2 + 2*r)*exp(-k*r/2)/(2*r*r_2d)
    
    psi = y*(-(*k)*z2 + 2*walker->get_r_i(i))/(2*walker->get_r_i(i)*(*r2d));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_3::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //y*z*(-4*r + r_2d^2*(k^2*r - 8*k))*exp(-k*r/2)/(4*r*r_2d^3)
    
    psi = y*z*(-4*walker->get_r_i(i) + (*r22d)*((*k2)*walker->get_r_i(i) - 8*(*k)))/(4*walker->get_r_i(i)*(*r22d)*(*r2d));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 3  -------------------------
*/

double hydrogenic_4::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //y*exp(-k*r/2)
    
    psi = y;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_4_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-k*x*y*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*x*y/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_4_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //(-k*y^2 + 2*r)*exp(-k*r/2)/(2*r)
    
    psi = (-(*k)*y2 + 2*walker->get_r_i(i))/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_4_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*y*z*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*y*z/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_4::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //k*y*(k*r - 8)*exp(-k*r/2)/(4*r)
    
    psi = (*k)*y*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 4  -------------------------
*/

