#ifndef GAUSSIANS_H
#define GAUSSIANS_H

#include "../Orbitals.h"

class Gaussians : public Orbitals
{
public:
    Gaussians(int n_p, const arma::uvec3 &powers, double alpha);

    // Orbitals interface
public:
    void set_qnum_indie_terms(Walker *walker, int i);
    double phi(const Walker *walker, int particle, int q_num);
    double del_phi(const Walker *walker, int particle, int q_num, int d);
    double lapl_phi(const Walker *walker, int particle, int q_num);

protected:
    double alpha;

    arma::uvec3 ABC;

    double* expFactor;
    arma::vec3 dFactor;

};

#endif // GAUSSIANS_H
