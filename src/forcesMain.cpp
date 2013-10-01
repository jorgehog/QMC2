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

void forceLoop(ASGD & asgd, VMC & vmc, DMC & dmc, SampleForce & forceSampler);

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
    minimizerParams.n_c_SGD=100;
    dmcParams.n_c=100;
    dmcParams.therm=200;
    dmcParams.n_b=20;
    generalParams.random_seed = -1380642306;

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

    if (parParams.is_master) std::cout << generalParams.random_seed << std::endl;

    //Creating a sampler
    SampleForce sampleForce(&generalParams.R, generalParams.n_p);

    //Creating the solver objects
    VMC vmc(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w);
    ASGD minimizer(&vmc, minimizerParams, parParams, generalParams.runpath);
    DMC dmc(generalParams, dmcParams, systemObjects, parParams);

    vmc.set_error_estimator(new SimpleVar(parParams));
    dmc.set_error_estimator(new SimpleVar(parParams));

    vmc.add_subsample(&sampleForce);
    dmc.add_subsample(&sampleForce);

    minimizer.minimize();
    vmc.run_method();

    double F = sampleForce.extract_mean();

    dmc.initFromVMC(&vmc);
    dmc.run_method();

    F = sampleForce.extract_mean_of_means();



    forceLoop(minimizer, vmc, dmc, sampleForce);


    if (parParams.is_master) cout << "~.* QMC fin *.~" << endl;

//Kinetic 1.221304
#ifdef MPI_ON
    MPI_Finalize();
#endif

    return 0;
}


void forceLoop(ASGD &asgd, VMC &vmc, DMC &dmc, SampleForce & forceSampler){

    using namespace arma;

    int nPoints = 20;
    double rMin = 0.5;
    double rMax = 5;

    vec rVec = linspace<vec>(rMin, rMax, nPoints);
    vec fVec = zeros<vec>(nPoints);

    for (int i = 0; i < nPoints; ++i) {

    }
}

#endif //QMC2_MAIN
