#include "nbodytransform.h"

#include "../../structs.h"
#include "../../Walker/Walker.h"
#include "../../Orbitals/OrbitalsFactory.h"

#include <iostream>


using namespace QMC2;


NBodyTransform::NBodyTransform(GeneralParams &gP, VariationalParams &vP,
                               const std::vector<BodyDef> &bodies,
                               OrbitalsFactory &factory) :
    Orbitals(gP.n_p, gP.dim),
    N(0)
{

    GeneralParams bodyParams;
    bodyParams.dim = gP.dim;

    int NTot = 0;

    populations.set_size(bodies.size());

    for (const BodyDef & body : bodies) {

        bodyParams.n_p = body.n_p_local;

        populations(N) = body.n_p_local;

        nuclei_walkers.push_back(new Walker(n_p, gP.dim));
        origins.push_back(body.origin);

        nuclei.push_back(factory.create(bodyParams, vP));
        NTot += body.n_p_local;
        N++;

    }

    name = TOSTR(N) + ("bodyTransformed" + nuclei.at(N-1)->getName());

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

    createMonoStructureCouplingMatrix();

    //Set the highest accessible state. Needed in order to calculate the
    //required exponentions (i.e. avoid calculating unneccessary ones)
    for (int i = 0; i < N; ++i) {
        nuclei.at(i)->set_nCap(2*(monoStructureCouplings.col(i).eval().max() + 1));
    }

}

void NBodyTransform::makePMatrix()
{
    P.set_size(n2, N);

    arma::mat tmp(N, N);
    tmp.ones();
    tmp.diag()*=-1;
    tmp(0, 0) = 1;

    int n_p_local;

    for (int i = 0; i < N; ++i) {
        n_p_local = populations(i);

        tmp.col(i) *= sqrt(n_p_local);
    }


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

void NBodyTransform::createMonoStructureCouplingMatrix()
{
    using namespace arma;

    monoStructureCouplings.set_size(n2, N);

#ifndef EXPERIMENTAL
    //
    for (int i = 0; i < n2; ++i) {
        for (int j = 0; j < N; ++j) {
            monoStructureCouplings(i,j) = i/N;
        }
    }
    return;
    //
#else

    unsigned int maxN2 = populations.max()/2;


    mat spEnergies(maxN2, N);
    for (unsigned int i = 0; i < maxN2; ++i) {
        for (int j = 0; j < N; ++j) {
            spEnergies(i, j) = nuclei.at(j)->get_sp_energy(i);
        }
    }

    int nCouplings = pow(maxN2, N);

    umat allQuantumNumbers(nCouplings, N);
    vec allEnergyCouplings(nCouplings);
    uvec indexMap(N);

    allEnergyCouplings.zeros();
    indexMap.zeros();

    for (int i = 0; i < nCouplings; ++i) {

        allQuantumNumbers.row(i) = indexMap.t();

        for (int j = 0; j < N; ++j) {
            allEnergyCouplings(i) += spEnergies(indexMap(j), j);
        }

        indexMap(N-1)++;

        for (int j = N-1; j >= 0; --j) {
            if (indexMap(j) == maxN2) {

                if (j == 0) {
                    break;
                }

                indexMap(j) = 0;
                indexMap(j-1)++;
            }
        }

    }

    uvec sortedIndexes = sort_index(allEnergyCouplings);

    int startRow;
    int endRow;
    double startE;
    double testE;

    uvec sumVec = sum(allQuantumNumbers, 1).eval()(sortedIndexes);
    uvec sumVecSpanSort;

    startRow = 0;
    endRow = 0;
    while (startRow < nCouplings - 1) {

        endRow = startRow + 1;

        startE = allEnergyCouplings(sortedIndexes(startRow));
        testE = allEnergyCouplings(sortedIndexes(endRow));

        if (startE != testE) {
            startRow = endRow;
            continue;
        }

        if (endRow < nCouplings - 1){
            while((allEnergyCouplings(sortedIndexes(endRow)) == startE) && (endRow < nCouplings - 1)) {
                endRow++;
            }
        }


        sumVecSpanSort = sort_index(sumVec(span(startRow, endRow)));
        sortedIndexes(span(startRow, endRow)) = sortedIndexes(span(startRow, endRow)).eval()(sumVecSpanSort);

        startRow = endRow;

    }

    for (int i = 0; i < n2; ++i) {
        monoStructureCouplings.row(i) = allQuantumNumbers.row(sortedIndexes(i/N));

    }
#endif
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


double NBodyTransform::get_dell_alpha_phi(Walker *walker, int p, int q_num, int n)
{
    (void) walker;

    double dell_alpha_phi = 0;

    for(int i = 0; i < N; i++){
        dell_alpha_phi += P(q_num, i)*nuclei.at(i)->get_dell_alpha_phi(nuclei_walkers.at(i), p, monoStructureCouplings(q_num, i), n);
    }

    return dell_alpha_phi;

}


double NBodyTransform::phi(const Walker *walker, int particle, int q_num)
{
    (void) walker;

    double phi = 0;

    for(int i = 0; i < N; i++){
        phi += P(q_num, i)*nuclei.at(i)->phi(nuclei_walkers.at(i), particle, monoStructureCouplings(q_num, i));
    }

    return phi;
}

double NBodyTransform::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    (void) walker;

    double del_phi = 0;

    for(int i = 0; i < N; i++){
        del_phi += P(q_num, i)*nuclei.at(i)->del_phi(nuclei_walkers.at(i), particle, monoStructureCouplings(q_num, i), d);
    }

    return del_phi;
}

double NBodyTransform::lapl_phi(const Walker *walker, int particle, int q_num)
{
    (void) walker;

    double lapl_phi = 0;

    for(int i = 0; i < N; i++){
        lapl_phi += P(q_num, i)*nuclei.at(i)->lapl_phi(nuclei_walkers.at(i), particle, monoStructureCouplings(q_num, i));
    }

    return lapl_phi;
}


void NBodyTransform::update(double R)
{
    double theta = 109.47; //deg
    origins.at(1) << R << 0 << 0;
    origins.at(2) << R*cos(theta) << R*sin(theta) << 0;
    makeRRelNucleiMatrix();
}

