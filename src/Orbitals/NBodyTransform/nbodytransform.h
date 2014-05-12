#pragma once


#include "../OrbitalsFactory.h"


namespace QMC2
{


struct GeneralParams;
struct VariationalParams;

struct BodyDef {
    int n_p_local;
    arma::rowvec origin;
};

class NBodyTransform : public Orbitals
{
public:
    NBodyTransform(GeneralParams & gP, VariationalParams &vP,
                   const std::vector<BodyDef> & bodies,
                   OrbitalsFactory &factory);

    friend class MolecularCoulomb;
protected:
    double get_parameter(int n);
    void set_parameter(double parameter, int n);
    double get_dell_alpha_phi(Walker *walker, int p, int q_num, int n);

public:
    void set_qnum_indie_terms(Walker *walker, int i);
    double phi(const Walker *walker, int particle, int q_num);
    double del_phi(const Walker *walker, int particle, int q_num, int d);
    double lapl_phi(const Walker *walker, int particle, int q_num);

    void debug();

    void update(double R);
private:

    int N;
    arma::mat P;

    std::vector<Walker*> nuclei_walkers;
    std::vector<Orbitals*> nuclei;
    arma::uvec populations;
    std::vector<arma::rowvec> origins;

    arma::mat r_rel_nuclei;
    arma::umat monoStructureCouplings;

    void makePMatrix();

    void makeRRelNucleiMatrix(); //-434.703048

    void createMonoStructureCouplingMatrix();

};

}
