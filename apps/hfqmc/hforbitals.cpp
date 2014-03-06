#include "hforbitals.h"


HFOrbitals::HFOrbitals(System *system) :
    m_system(system)
{

}

double HFOrbitals::phi(const QMC2::Walker *walker, int particle, int q_num)
{
    const rowvec & ri = walker->r.row(particle);

    return m_system->evaluateCGTO(q_num, ri(0), ri(1), ri(2));
}

double HFOrbitals::del_phi(const QMC2::Walker *walker, int particle, int q_num, int d)
{
    num_diff(walker, particle, q_num, d);
}

double HFOrbitals::lapl_phi(const QMC2::Walker *walker, int particle, int q_num)
{
    num_ddiff(walker, particle, q_num);
}
