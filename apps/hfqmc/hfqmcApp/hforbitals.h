#pragma once

#include <QMC2.h>
#include <hf.h>

namespace QMC2
{

class HFOrbitals : public Orbitals
{
public:
    HFOrbitals(const uint n_p, vector<const hf::ContractedGTO *> contractedGTOs, const mat &coeffs);

    double evalContracted(uint k, const rowvec &r);
    double evallaplContracted(uint k, const rowvec &r);
    double evaldelContracted(uint k, const rowvec &r, int d);

    // Orbitals interface
public:
    void set_qnum_indie_terms(Walker *walker, int i);
    double phi(const Walker *walker, int particle, int q_num);
    double del_phi(const Walker *walker, int particle, int q_num, int d);
    double lapl_phi(const Walker *walker, int particle, int q_num);

    vector<const hf::ContractedGTO *>  m_contractedGTOs;
    const mat m_coeffs;


};

}
