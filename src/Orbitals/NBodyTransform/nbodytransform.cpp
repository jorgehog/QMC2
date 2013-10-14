#include "nbodytransform.h"

#include "../../structs.h"

#include "../../Walker/Walker.h"

#include "../AlphaHarmonicOscillator/AlphaHarmonicOscillator.h"
#include "../hydrogenicOrbitals/hydrogenicOrbitals.h"

#include <iostream>

NBodyTransform::NBodyTransform(GeneralParams &gP, VariationalParams &vP, TRANS_SYSTEMS system, const std::vector<BodyDef> &bodies) :
    Orbitals(gP.n_p, gP.dim),
    N(0)
{

    GeneralParams bodyParams;
    bodyParams.dim = gP.dim;

    int NTot = 0;

    for (const BodyDef & body : bodies) {
        bodyParams.n_p = body.n_p_local;

        populations.push_back(body.n_p_local);

        nuclei_walkers.push_back(new Walker(n_p, gP.dim));
        origins.push_back(body.origin);

        switch (system) {
        case ATOMS:

            name = "Molecule";

            nuclei.push_back(new hydrogenicOrbitals(bodyParams, vP));
            break;

        case QDOTS:

            name = "QDotLattice";

            nuclei.push_back(new AlphaHarmonicOscillator(bodyParams, vP));
            break;

        default:

            std::cout << "System valued" << system << " has not been implemented for DiTransform." << std::endl;
            exit(1);

            break;
        }

        NTot += body.n_p_local;
        N++;

    }

    if (N <= 1) {
        std::cout << "No bodies transformed..." << std::endl;
        exit(1);
    }
    if (NTot != n_p) {
        std::cout << "nBodyTransform failed. Mismatch in the particle numbers." << std::endl;
        exit(1);
    }

    makeRRelNucleiMatrix();

    makePMatrix();

}

void NBodyTransform::makePMatrix()
{
    P.set_size(n2, N);

    arma::mat tmp(N, N);
    tmp.ones();
    tmp.diag()*=-1;
    tmp(0, 0) = 1;

//    int nFac = 0;
    int n_p_local;

    for (int i = 0; i < N; ++i) {
        n_p_local = populations.at(i);

        tmp.col(i) *= ((double)n_p_local)/n_p;

//        nFac += n_p_local*n_p_local;
    }

//    tmp *= n_p/sqrt(nFac);

    for (int i = 0; i < n2; ++i) {

        P.row(i) = tmp.row(i%N);

    }

}

void NBodyTransform::makeRRelNucleiMatrix()
{
    r_rel_nuclei.set_size(N, N);

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            r_rel_nuclei(i, j) = r_rel_nuclei(j, i) = arma::norm(origins.at(i) - origins.at(j), 2);
        }
    }

}

void NBodyTransform::set_qnum_indie_terms(Walker *walker, int i)
{
//    walker->calc_r_i(i);

    for (int k = 0; k < N; ++k) {

        nuclei_walkers.at(k)->r.row(i) = walker->r.row(i) - origins.at(k);

        nuclei_walkers.at(k)->calc_r_i2(i);

        nuclei.at(k)->set_qnum_indie_terms(nuclei_walkers.at(k), i);

    }

}

double NBodyTransform::get_parameter(int n)
{
    return (*(nuclei.begin()))->get_parameter(n);
}


void NBodyTransform::set_parameter(double parameter, int n)
{
    for (Orbitals* nucleus : nuclei) {
        nucleus->set_parameter(parameter, n);
    }
}


double NBodyTransform::get_dell_alpha_phi(Walker *walker, int p, int q_num)
{
    (void) walker;

    double dell_alpha_phi = 0;

    for(int i = 0; i < N; i++){
        dell_alpha_phi += P(q_num, i)*nuclei.at(i)->get_dell_alpha_phi(nuclei_walkers.at(i), p, q_num/N);
    }

    return dell_alpha_phi;

}


double NBodyTransform::phi(const Walker *walker, int particle, int q_num)
{
    (void) walker;

    double phi = 0;

    for(int i = 0; i < N; i++){
        phi += P(q_num, i)*nuclei.at(i)->phi(nuclei_walkers.at(i), particle, q_num/N);
    }

    return phi;
}

double NBodyTransform::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    (void) walker;

    double del_phi = 0;

    for(int i = 0; i < N; i++){
        del_phi += P(q_num, i)*nuclei.at(i)->del_phi(nuclei_walkers.at(i), particle, q_num/N, d);
    }

    return del_phi;
}

double NBodyTransform::lapl_phi(const Walker *walker, int particle, int q_num)
{
    (void) walker;

    double lapl_phi = 0;

    for(int i = 0; i < N; i++){
        lapl_phi += P(q_num, i)*nuclei.at(i)->lapl_phi(nuclei_walkers.at(i), particle, q_num/N);
    }

    return lapl_phi;
}

void NBodyTransform::debug()
{
    std::cout << *((hydrogenicOrbitals*)nuclei.at(0))->exp_factor_n1 << std::endl;
    std::cout << *((hydrogenicOrbitals*)nuclei.at(1))->exp_factor_n1 << std::endl;
    std::cout << *((hydrogenicOrbitals*)nuclei.at(1))->k << std::endl;
}

