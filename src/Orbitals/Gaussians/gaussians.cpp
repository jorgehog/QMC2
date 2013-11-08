#include "gaussians.h"

#include "../../Walker/Walker.h"

Gaussians::Gaussians(int n_p, const arma::uvec3 &powers, double alpha) :
    alpha(alpha),
    ABC(powers)
{
    set_n_p(n_p);
    set_dim(3);
}


void Gaussians::set_qnum_indie_terms(Walker *walker, int i)
{

    *expFactor = exp(-alpha*walker->get_r_i2(i));

    for (int d = 0; d < 3; ++d) {
        const double & xi = walker->r(i, d);
        dFactor(d) = (ABC(d)/xi - 2*alpha*xi);
    }

}

double Gaussians::phi(const Walker *walker, int particle, int q_num)
{

    (void) q_num;

    double G = 1;
    for (int i = 0; i < 3; ++i) {
        G *= pow(walker->r(particle, i), ABC(i));
    }
    return G*(*expFactor);
}

double Gaussians::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    return walker->phi(particle, q_num)*dFactor(d);
}

double Gaussians::lapl_phi(const Walker *walker, int particle, int q_num)
{

    (void) q_num;

    double ddG = 0;

    for (int d = 0; d < 3; ++d) {

        const double & xi = walker->r(particle, d);

        ddG += dFactor(d)*dFactor(d) - ABC(d)/(xi*xi);

    }

    ddG -= 6*alpha;

    return walker->phi(particle, q_num)*ddG;

}
