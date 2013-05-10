
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

hydrogenic_5::hydrogenic_5(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_5_x::dell_hydrogenic_5_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_5_y::dell_hydrogenic_5_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_5_z::dell_hydrogenic_5_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_5::lapl_hydrogenic_5(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 5  -------------------------
*/

hydrogenic_6::hydrogenic_6(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_6_x::dell_hydrogenic_6_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_6_y::dell_hydrogenic_6_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_6_z::dell_hydrogenic_6_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_6::lapl_hydrogenic_6(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 6  -------------------------
*/

hydrogenic_7::hydrogenic_7(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_7_x::dell_hydrogenic_7_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_7_y::dell_hydrogenic_7_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_7_z::dell_hydrogenic_7_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_7::lapl_hydrogenic_7(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 7  -------------------------
*/

hydrogenic_8::hydrogenic_8(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_8_x::dell_hydrogenic_8_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_8_y::dell_hydrogenic_8_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_8_z::dell_hydrogenic_8_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_8::lapl_hydrogenic_8(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 8  -------------------------
*/

hydrogenic_9::hydrogenic_9(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_9_x::dell_hydrogenic_9_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_9_y::dell_hydrogenic_9_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_9_z::dell_hydrogenic_9_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_9::lapl_hydrogenic_9(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 9  -------------------------
*/

hydrogenic_10::hydrogenic_10(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_10_x::dell_hydrogenic_10_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_10_y::dell_hydrogenic_10_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_10_z::dell_hydrogenic_10_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_10::lapl_hydrogenic_10(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 10  -------------------------
*/

hydrogenic_11::hydrogenic_11(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_11_x::dell_hydrogenic_11_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_11_y::dell_hydrogenic_11_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_11_z::dell_hydrogenic_11_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_11::lapl_hydrogenic_11(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 11  -------------------------
*/

hydrogenic_12::hydrogenic_12(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_12_x::dell_hydrogenic_12_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_12_y::dell_hydrogenic_12_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_12_z::dell_hydrogenic_12_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_12::lapl_hydrogenic_12(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 12  -------------------------
*/

hydrogenic_13::hydrogenic_13(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_13_x::dell_hydrogenic_13_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_13_y::dell_hydrogenic_13_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_13_z::dell_hydrogenic_13_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_13::lapl_hydrogenic_13(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 13  -------------------------
*/

hydrogenic_14::hydrogenic_14(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_14_x::dell_hydrogenic_14_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_14_y::dell_hydrogenic_14_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_14_z::dell_hydrogenic_14_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_14::lapl_hydrogenic_14(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 14  -------------------------
*/

hydrogenic_15::hydrogenic_15(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_15_x::dell_hydrogenic_15_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_15_y::dell_hydrogenic_15_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_15_z::dell_hydrogenic_15_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_15::lapl_hydrogenic_15(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 15  -------------------------
*/

hydrogenic_16::hydrogenic_16(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_16_x::dell_hydrogenic_16_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_16_y::dell_hydrogenic_16_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_16_z::dell_hydrogenic_16_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_16::lapl_hydrogenic_16(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 16  -------------------------
*/

hydrogenic_17::hydrogenic_17(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_17_x::dell_hydrogenic_17_x(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_17_y::dell_hydrogenic_17_y(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

dell_hydrogenic_17_z::dell_hydrogenic_17_z(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

lapl_hydrogenic_17::lapl_hydrogenic_17(double* k, double* k2, double* r22d, double* r2d, double* exp_factor)
: hydrogenic(k, k2, r22d, r2d, exp_factor) {

}

/*
    -------------------------  END 17  -------------------------
*/



/*
    Subclass Eval functions
*/


double hydrogenic_0::eval(const Walker* walker, int i) {
    
    //exp(-k*r)
    
    psi = 1;
    //std::cout << "hydro 0  " << psi << " " << *exp_factor << std::endl;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_0_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //-k*x*exp(-k*r)/r
    
    psi = -(*k)*x/walker->get_r_i(i);
    //std::cout << "hydro 0_dx  " << psi << " " << *exp_factor << std::endl;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_0_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //-k*y*exp(-k*r)/r
    
    psi = -(*k)*y/walker->get_r_i(i);
    //std::cout << "hydro 0_dy  " << psi << " " << *exp_factor << std::endl;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_0_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //-k*z*exp(-k*r)/r
    
    psi = -(*k)*z/walker->get_r_i(i);
    //std::cout << "hydro 0_dz  " << psi << " " << *exp_factor << std::endl;
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_0::eval(const Walker* walker, int i) {
    
    //k*(k*r - 2)*exp(-k*r)/r
    
    psi = (*k)*((*k)*walker->get_r_i(i) - 2)/walker->get_r_i(i);
    //std::cout << "hydro 0_lapl  " << psi << " " << *exp_factor << std::endl;
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

    z = walker->r(i, 2);
    
    //z*exp(-k*r/2)
    
    psi = z;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_2_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-k*x*z*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*x*z/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_2_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*y*z*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*y*z/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_2_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //(-k*z^2 + 2*r)*exp(-k*r/2)/(2*r)
    
    psi = (-(*k)*z2 + 2*walker->get_r_i(i))/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_2::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //k*z*(k*r - 8)*exp(-k*r/2)/(4*r)
    
    psi = (*k)*z*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 2  -------------------------
*/

double hydrogenic_3::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //x*exp(-k*r/2)
    
    psi = x;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_3_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //(-k*x^2 + 2*r)*exp(-k*r/2)/(2*r)
    
    psi = (-(*k)*x2 + 2*walker->get_r_i(i))/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_3_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-k*x*y*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*x*y/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_3_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-k*x*z*exp(-k*r/2)/(2*r)
    
    psi = -(*k)*x*z/(2*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_3::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //k*x*(k*r - 8)*exp(-k*r/2)/(4*r)
    
    psi = (*k)*x*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
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

double hydrogenic_5::eval(const Walker* walker, int i) {
    
    //(2*k^2*r^2 - 18*k*r + 27)*exp(-k*r/3)
    
    psi = 2*(*k2)*walker->get_r_i2(i) - 18*(*k)*walker->get_r_i(i) + 27;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_5_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //-k*x*(2*k^2*r^2 - 30*k*r + 81)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*(2*(*k2)*walker->get_r_i2(i) - 30*(*k)*walker->get_r_i(i) + 81)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_5_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //-k*y*(2*k^2*r^2 - 30*k*r + 81)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*y*(2*(*k2)*walker->get_r_i2(i) - 30*(*k)*walker->get_r_i(i) + 81)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_5_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //-k*z*(2*k^2*r^2 - 30*k*r + 81)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*z*(2*(*k2)*walker->get_r_i2(i) - 30*(*k)*walker->get_r_i(i) + 81)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_5::eval(const Walker* walker, int i) {
    
    //k*(k*r - 18)*(2*k^2*r^2 - 18*k*r + 27)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*((*k)*walker->get_r_i(i) - 18)*(2*(*k2)*walker->get_r_i2(i) - 18*(*k)*walker->get_r_i(i) + 27)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 5  -------------------------
*/

double hydrogenic_6::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //z*(k*r - 6)*exp(-k*r/3)
    
    psi = z*((*k)*walker->get_r_i(i) - 6);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_6_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-k*x*z*(k*r - 9)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*z*((*k)*walker->get_r_i(i) - 9)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_6_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*y*z*(k*r - 9)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*y*z*((*k)*walker->get_r_i(i) - 9)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_6_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //(3*k*r^2 - k*z^2*(k*r - 9) - 18*r)*exp(-k*r/3)/(3*r)
    
    psi = (3*(*k)*walker->get_r_i2(i) - (*k)*z2*((*k)*walker->get_r_i(i) - 9) - 18*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_6::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //k*z*(k*r - 18)*(k*r - 6)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*z*((*k)*walker->get_r_i(i) - 18)*((*k)*walker->get_r_i(i) - 6)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 6  -------------------------
*/

double hydrogenic_7::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //x*(k*r - 6)*exp(-k*r/3)
    
    psi = x*((*k)*walker->get_r_i(i) - 6);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_7_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //(3*k*r^2 - k*x^2*(k*r - 9) - 18*r)*exp(-k*r/3)/(3*r)
    
    psi = (3*(*k)*walker->get_r_i2(i) - (*k)*x2*((*k)*walker->get_r_i(i) - 9) - 18*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_7_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-k*x*y*(k*r - 9)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*y*((*k)*walker->get_r_i(i) - 9)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_7_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-k*x*z*(k*r - 9)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*z*((*k)*walker->get_r_i(i) - 9)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_7::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //k*x*(k*r - 18)*(k*r - 6)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*x*((*k)*walker->get_r_i(i) - 18)*((*k)*walker->get_r_i(i) - 6)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 7  -------------------------
*/

double hydrogenic_8::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //y*(k*r - 6)*exp(-k*r/3)
    
    psi = y*((*k)*walker->get_r_i(i) - 6);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_8_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-k*x*y*(k*r - 9)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*y*((*k)*walker->get_r_i(i) - 9)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_8_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //(3*k*r^2 - k*y^2*(k*r - 9) - 18*r)*exp(-k*r/3)/(3*r)
    
    psi = (3*(*k)*walker->get_r_i2(i) - (*k)*y2*((*k)*walker->get_r_i(i) - 9) - 18*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_8_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*y*z*(k*r - 9)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*y*z*((*k)*walker->get_r_i(i) - 9)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_8::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //k*y*(k*r - 18)*(k*r - 6)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*y*((*k)*walker->get_r_i(i) - 18)*((*k)*walker->get_r_i(i) - 6)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 8  -------------------------
*/

double hydrogenic_9::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //(-r^2 + 3*z^2)*exp(-k*r/3)
    
    psi = -walker->get_r_i2(i) + 3*z2;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_9_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-x*(k*(-r^2 + 3*z^2) + 6*r)*exp(-k*r/3)/(3*r)
    
    psi = -x*((*k)*(-walker->get_r_i2(i) + 3*z2) + 6*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_9_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-y*(k*(-r^2 + 3*z^2) + 6*r)*exp(-k*r/3)/(3*r)
    
    psi = -y*((*k)*(-walker->get_r_i2(i) + 3*z2) + 6*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_9_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //-z*(k*(-r^2 + 3*z^2) - 12*r)*exp(-k*r/3)/(3*r)
    
    psi = -z*((*k)*(-walker->get_r_i2(i) + 3*z2) - 12*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_9::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //k*(-r^2 + 3*z^2)*(k*r - 18)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*(-walker->get_r_i2(i) + 3*z2)*((*k)*walker->get_r_i(i) - 18)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 9  -------------------------
*/

double hydrogenic_10::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //x*z*exp(-k*r/3)
    
    psi = x*z;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_10_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    x2 = x*x;
    
    //-z*(k*x^2 - 3*r)*exp(-k*r/3)/(3*r)
    
    psi = -z*((*k)*x2 - 3*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_10_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*x*y*z*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*y*z/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_10_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-x*(k*z^2 - 3*r)*exp(-k*r/3)/(3*r)
    
    psi = -x*((*k)*z2 - 3*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_10::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //k*x*z*(k*r - 18)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*x*z*((*k)*walker->get_r_i(i) - 18)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 10  -------------------------
*/

double hydrogenic_11::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //y*z*exp(-k*r/3)
    
    psi = y*z;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_11_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*x*y*z*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*y*z/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_11_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    y2 = y*y;
    
    //-z*(k*y^2 - 3*r)*exp(-k*r/3)/(3*r)
    
    psi = -z*((*k)*y2 - 3*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_11_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);

    z2 = z*z;
    
    //-y*(k*z^2 - 3*r)*exp(-k*r/3)/(3*r)
    
    psi = -y*((*k)*z2 - 3*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_11::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //k*y*z*(k*r - 18)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*y*z*((*k)*walker->get_r_i(i) - 18)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 11  -------------------------
*/

double hydrogenic_12::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //(x^2 - y^2)*exp(-k*r/3)
    
    psi = x2 - y2;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_12_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //-x*(k*(x^2 - y^2) - 6*r)*exp(-k*r/3)/(3*r)
    
    psi = -x*((*k)*(x2 - y2) - 6*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_12_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //-y*(k*(x^2 - y^2) + 6*r)*exp(-k*r/3)/(3*r)
    
    psi = -y*((*k)*(x2 - y2) + 6*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_12_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);

    x2 = x*x;
    y2 = y*y;
    
    //-k*z*(x^2 - y^2)*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*z*(x2 - y2)/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_12::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    y2 = y*y;
    
    //k*(x^2 - y^2)*(k*r - 18)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*(x2 - y2)*((*k)*walker->get_r_i(i) - 18)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 12  -------------------------
*/

double hydrogenic_13::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //x*y*exp(-k*r/3)
    
    psi = x*y;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_13_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    x2 = x*x;
    
    //-y*(k*x^2 - 3*r)*exp(-k*r/3)/(3*r)
    
    psi = -y*((*k)*x2 - 3*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_13_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);

    y2 = y*y;
    
    //-x*(k*y^2 - 3*r)*exp(-k*r/3)/(3*r)
    
    psi = -x*((*k)*y2 - 3*walker->get_r_i(i))/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_13_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*x*y*z*exp(-k*r/3)/(3*r)
    
    psi = -(*k)*x*y*z/(3*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_13::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //k*x*y*(k*r - 18)*exp(-k*r/3)/(9*r)
    
    psi = (*k)*x*y*((*k)*walker->get_r_i(i) - 18)/(9*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 13  -------------------------
*/

double hydrogenic_14::eval(const Walker* walker, int i) {
    
    //(k^3*r^3 - 24*k^2*r^2 + 144*k*r - 192)*exp(-k*r/4)
    
    psi = (*k2)*(*k)*walker->get_r_i2(i)*walker->get_r_i(i) - 24*(*k2)*walker->get_r_i2(i) + 144*(*k)*walker->get_r_i(i) - 192;
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_14_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //-k*x*(k^3*r^3 - 36*k^2*r^2 + 336*k*r - 768)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*x*((*k2)*(*k)*walker->get_r_i2(i)*walker->get_r_i(i) - 36*(*k2)*walker->get_r_i2(i) + 336*(*k)*walker->get_r_i(i) - 768)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_14_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //-k*y*(k^3*r^3 - 36*k^2*r^2 + 336*k*r - 768)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*y*((*k2)*(*k)*walker->get_r_i2(i)*walker->get_r_i(i) - 36*(*k2)*walker->get_r_i2(i) + 336*(*k)*walker->get_r_i(i) - 768)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_14_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //-k*z*(k^3*r^3 - 36*k^2*r^2 + 336*k*r - 768)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*z*((*k2)*(*k)*walker->get_r_i2(i)*walker->get_r_i(i) - 36*(*k2)*walker->get_r_i2(i) + 336*(*k)*walker->get_r_i(i) - 768)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_14::eval(const Walker* walker, int i) {
    
    //k*(k*r - 32)*(k^3*r^3 - 24*k^2*r^2 + 144*k*r - 192)*exp(-k*r/4)/(16*r)
    
    psi = (*k)*((*k)*walker->get_r_i(i) - 32)*((*k2)*(*k)*walker->get_r_i2(i)*walker->get_r_i(i) - 24*(*k2)*walker->get_r_i2(i) + 144*(*k)*walker->get_r_i(i) - 192)/(16*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 14  -------------------------
*/

double hydrogenic_15::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //z*(k^2*r^2 - 20*k*r + 80)*exp(-k*r/4)
    
    psi = z*((*k2)*walker->get_r_i2(i) - 20*(*k)*walker->get_r_i(i) + 80);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_15_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-k*x*z*(k*r - 20)*(k*r - 8)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*x*z*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_15_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*y*z*(k*r - 20)*(k*r - 8)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*y*z*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_15_z::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);

    z2 = z*z;
    
    //(4*k^2*r^3 - 80*k*r^2 - k*z^2*(k*r - 20)*(k*r - 8) + 320*r)*exp(-k*r/4)/(4*r)
    
    psi = (4*(*k2)*walker->get_r_i2(i)*walker->get_r_i(i) - 80*(*k)*walker->get_r_i2(i) - (*k)*z2*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8) + 320*walker->get_r_i(i))/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_15::eval(const Walker* walker, int i) {

    z = walker->r(i, 2);
    
    //k*z*(k*r - 32)*(k^2*r^2 - 20*k*r + 80)*exp(-k*r/4)/(16*r)
    
    psi = (*k)*z*((*k)*walker->get_r_i(i) - 32)*((*k2)*walker->get_r_i2(i) - 20*(*k)*walker->get_r_i(i) + 80)/(16*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 15  -------------------------
*/

double hydrogenic_16::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //x*(k^2*r^2 - 20*k*r + 80)*exp(-k*r/4)
    
    psi = x*((*k2)*walker->get_r_i2(i) - 20*(*k)*walker->get_r_i(i) + 80);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_16_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);

    x2 = x*x;
    
    //(4*k^2*r^3 - 80*k*r^2 - k*x^2*(k*r - 20)*(k*r - 8) + 320*r)*exp(-k*r/4)/(4*r)
    
    psi = (4*(*k2)*walker->get_r_i2(i)*walker->get_r_i(i) - 80*(*k)*walker->get_r_i2(i) - (*k)*x2*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8) + 320*walker->get_r_i(i))/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_16_y::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-k*x*y*(k*r - 20)*(k*r - 8)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*x*y*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_16_z::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    z = walker->r(i, 2);
    
    //-k*x*z*(k*r - 20)*(k*r - 8)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*x*z*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_16::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    
    //k*x*(k*r - 32)*(k^2*r^2 - 20*k*r + 80)*exp(-k*r/4)/(16*r)
    
    psi = (*k)*x*((*k)*walker->get_r_i(i) - 32)*((*k2)*walker->get_r_i2(i) - 20*(*k)*walker->get_r_i(i) + 80)/(16*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 16  -------------------------
*/

double hydrogenic_17::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //y*(k^2*r^2 - 20*k*r + 80)*exp(-k*r/4)
    
    psi = y*((*k2)*walker->get_r_i2(i) - 20*(*k)*walker->get_r_i(i) + 80);
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_17_x::eval(const Walker* walker, int i) {

    x = walker->r(i, 0);
    y = walker->r(i, 1);
    
    //-k*x*y*(k*r - 20)*(k*r - 8)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*x*y*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_17_y::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);

    y2 = y*y;
    
    //(4*k^2*r^3 - 80*k*r^2 - k*y^2*(k*r - 20)*(k*r - 8) + 320*r)*exp(-k*r/4)/(4*r)
    
    psi = (4*(*k2)*walker->get_r_i2(i)*walker->get_r_i(i) - 80*(*k)*walker->get_r_i2(i) - (*k)*y2*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8) + 320*walker->get_r_i(i))/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double dell_hydrogenic_17_z::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    z = walker->r(i, 2);
    
    //-k*y*z*(k*r - 20)*(k*r - 8)*exp(-k*r/4)/(4*r)
    
    psi = -(*k)*y*z*((*k)*walker->get_r_i(i) - 20)*((*k)*walker->get_r_i(i) - 8)/(4*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

double lapl_hydrogenic_17::eval(const Walker* walker, int i) {

    y = walker->r(i, 1);
    
    //k*y*(k*r - 32)*(k^2*r^2 - 20*k*r + 80)*exp(-k*r/4)/(16*r)
    
    psi = (*k)*y*((*k)*walker->get_r_i(i) - 32)*((*k2)*walker->get_r_i2(i) - 20*(*k)*walker->get_r_i(i) + 80)/(16*walker->get_r_i(i));
    return psi*(*exp_factor);
    
}

/*
    -------------------------  END 17  -------------------------
*/

