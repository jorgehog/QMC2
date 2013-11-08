#ifndef GASSIANSET_H
#define GASSIANSET_H

#include "../Orbitals.h"

class Gaussians;

class GassianSet : public Orbitals
{
public:
    GassianSet(int n_p, arma::vec D, arma::uvec3 powers, arma::vec3 alpha);

    // Orbitals interface
public:
    void set_qnum_indie_terms(Walker *walker, int i);
    double phi(const Walker *walker, int particle, int q_num);
    double del_phi(const Walker *walker, int particle, int q_num, int d);
    double lapl_phi(const Walker *walker, int particle, int q_num);
private:

    std::vector<Gaussians*> gaussians;
    arma::vec D;

};

#endif // GASSIANSET_H
