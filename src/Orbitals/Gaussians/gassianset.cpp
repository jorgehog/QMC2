#include "gassianset.h"

#include "gaussians.h"

GassianSet::GassianSet(int n_p, arma::vec D, arma::uvec3 powers, arma::vec3 alpha) :
    D(D)
{

    set_n_p(n_p);
    set_dim(3);

    for (int i = 0; i < D.n_elem; ++i) {
        gaussians.push_back(new Gaussians(n_p, powers, alpha(i)));
    }

}


void GassianSet::set_qnum_indie_terms(Walker *walker, int i)
{
    for(Gaussians * gaussian : gaussians) {
        gaussian->set_qnum_indie_terms(walker, i);
    }

}

double GassianSet::phi(const Walker *walker, int particle, int q_num)
{
    double _phi = 0;
    for (int i = 0; i < gaussians.size(); ++i) {
        _phi += D(i)*gaussians.at(i)->phi(walker, particle, q_num);
    }

    return _phi;
}

double GassianSet::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    double _phi = 0;
    for (int i = 0; i < gaussians.size(); ++i) {
        _phi += D(i)*gaussians.at(i)->del_phi(walker, particle, q_num, d);
    }

    return _phi;
}

double GassianSet::lapl_phi(const Walker *walker, int particle, int q_num)
{
    double _phi = 0;
    for (int i = 0; i < gaussians.size(); ++i) {
        _phi += D(i)*gaussians.at(i)->lapl_phi(walker, particle, q_num);
    }

    return _phi;
}
