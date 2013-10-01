/*
 * File:   QMCmain.cpp
 * Author: jorgehog
 *
 * Created on 13. april 2012, 17:04
 */

#ifndef QMC2_MAIN
#define QMC2_MAIN

//#include "HartreeFock/HartreeFock.h"

#include "structs.h"

#include "QMC/VMC/VMC.h"
#include "QMC/DMC/DMC.h"

#include "Minimizer/ASGD/ASGD.h"

#include "ErrorEstimator/SimpleVar/SimpleVar.h"
#include "ErrorEstimator/Blocking/Blocking.h"

#include "OutputHandler/stdoutASGD/stdoutASGD.h"
#include "OutputHandler/stdoutDMC/stdoutDMC.h"
#include "OutputHandler/Distribution/Distribution.h"

#include "Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h"
#include "Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.h"
#include "Orbitals/DiTransform/DiTransform.h"
//#include "Orbitals/ExpandedBasis/ExpandedBasis.h"

#include "Potential/AtomCore/AtomCore.h"
#include "Potential/Coulomb/Coulomb.h"
#include "Potential/DiAtomCore/DiAtomCore.h"
#include "Potential/DoubleWell/DoubleWell.h"
#include "Potential/Harmonic_osc/Harmonic_osc.h"

//#include "System/Bosons/Bosons.h"
#include "System/Fermions/Fermions.h"

#include "Sampling/Brute_Force/Brute_Force.h"
#include "Sampling/Importance/Importance.h"

#include "Jastrow/No_Jastrow/No_Jastrow.h"
#include "Jastrow/Pade_Jastrow/Pade_Jastrow.h"

#include "Sampler/sampleMethods/SampleForce.h"

void forceLoop(ASGD & asgd, VMC & vmc, DMC & dmc, double *R, int n_p, bool is_master);

int main(int argc, char** argv) {

    using namespace std;

    //Initialize structs for holding constructor arguments.
    struct ParParams parParams;
    struct VMCparams vmcParams;
    struct DMCparams dmcParams;
    struct VariationalParams variationalParams;
    struct GeneralParams generalParams;
    struct MinimizerParams minimizerParams;
    struct SystemObjects systemObjects;

    minimizerParams.n_c_SGD=600;
    minimizerParams.SGDsamples = 4000;
    dmcParams.n_c=250;
    dmcParams.therm=1000;
    dmcParams.n_b=50;
    dmcParams.n_w = 2000;

    //Setting up parallel parameters
    initMPI(parParams, argc, argv);
    scaleWithProcs(parParams, generalParams, minimizerParams, vmcParams, dmcParams);

    //Setting up the solver parameters
    generalParams.dim = 3;
    systemObjects.SP_basis = new DiTransform(generalParams, variationalParams, ATOMS);

    Fermions system(generalParams, systemObjects.SP_basis);
    system.add_potential(new DiAtomCore(generalParams));
    system.add_potential(new Coulomb(generalParams));
    systemObjects.system = &system;

    systemObjects.jastrow = new Pade_Jastrow(generalParams, variationalParams);
    systemObjects.sample_method = new Importance(generalParams);

    //

    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;

    //Creating the solver objects
    bool silent = true;
    VMC vmc(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w, silent);
    ASGD minimizer(&vmc, minimizerParams, parParams, generalParams.runpath);
    DMC dmc(generalParams, dmcParams, systemObjects, parParams, silent);

    vmc.set_error_estimator(new SimpleVar(parParams));
    dmc.set_error_estimator(new SimpleVar(parParams));

    forceLoop(minimizer, vmc, dmc, &generalParams.R, generalParams.n_p, parParams.is_master);

    if (parParams.is_master) cout << "~.* QMC fin *.~" << endl;

#ifdef MPI_ON
    MPI_Finalize();
#endif

    return 0;
}


void forceLoop(ASGD &asgd, VMC &vmc, DMC &dmc, double *R, int n_p, bool is_master){

    using namespace arma;

    double f_dmc, f_vmc, f0;

    int nPoints = 100;
    double rMin = 0.5;
    double rMax = 5;

    vec rVec = linspace<vec>(rMin, rMax, nPoints);
    vec fVec = zeros<vec>(nPoints);

    SampleForce forceSampler(R, n_p);
    vmc.add_subsample(&forceSampler);
    dmc.add_subsample(&forceSampler);



    //Performing first iteration with initializations
    *R = rMin;

    asgd.minimize();
    vmc.run_method();

    f_vmc = forceSampler.extract_mean();

    dmc.initFromVMC(&vmc);
    dmc.run_method();

    f_dmc = forceSampler.extract_mean_of_means();

    f0 = 2*f_dmc - f_vmc;
    fVec(0) = f0;

    fVec.save("/tmp/binaryArmaVec.arma");

    for (int i = 1; i < nPoints; ++i) {

        *R = rVec(i);

        if (is_master) std::cout << "\nLooping R = " << *R << "\n" << std::endl;

        asgd.minimize(false);
        vmc.run_method(false);

        f_vmc = forceSampler.extract_mean();

        dmc.run_method(false);

        f_dmc = forceSampler.extract_mean_of_means();

        f0 = 2*f_dmc - f_vmc;
        fVec(i) = f0;

        fVec.save("/tmp/binaryArmaVec.arma");

    }
}

#endif //QMC2_MAIN
