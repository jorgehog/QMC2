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

#include "Orbitals/OrbitalsFactory.h"
#include "Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h"
#include "Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.h"
#include "Orbitals/Gaussians/oxygen3-21G/oxygen3_21g.h"
#include "Orbitals/NBodyTransform/nbodytransform.h"
#include "Orbitals/ExpandedBasis/ExpandedBasis.h"

#include "Potential/AtomCore/AtomCore.h"
#include "Potential/Coulomb/Coulomb.h"
#include "Potential/DoubleWell/DoubleWell.h"
#include "Potential/Harmonic_osc/Harmonic_osc.h"
#include "Potential/MolecularCoulomb/molecularcoulomb.h"

//#include "System/Bosons/Bosons.h"
#include "System/Fermions/Fermions.h"

#include "Sampling/Brute_Force/Brute_Force.h"
#include "Sampling/Importance/Importance.h"

#include "Jastrow/No_Jastrow/No_Jastrow.h"
#include "Jastrow/Pade_Jastrow/Pade_Jastrow.h"

#include "Sampler/sampleMethods/SampleForce.h"

#include <armadillo>
#include <iomanip>
#include <vector>

void forceLoop(ASGD & asgd, VMC & vmc, DMC & dmc, int n_p, bool is_master);

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

    if (argc == 2) {
        generalParams.runpath = argv[1];
    }

    minimizerParams.SGDsamples = 3000;

    dmcParams.therm=1000;
    dmcParams.n_c = 2000;
    vmcParams.n_c = 1E5;
    dmcParams.dt = 0.00005;
    generalParams.random_seed = -1384354898;

    variationalParams.alpha = 0.922925;
    variationalParams.beta  = 0.3477;

    //Setting up parallel parameters
    initMPI(parParams, argc, argv);

    scaleWithProcs(parParams, generalParams, minimizerParams, vmcParams, dmcParams);

    //Setting up the solver parameters
    generalParams.n_p = 2;
    generalParams.dim = 3;

//    gP.dim = 3;

//    sO.SP_basis = new hydrogenicOrbitals(gP, vP);

//    sO.onebody_pot = new AtomCore(gP);

//    sO.system = new Fermions(gP, sO.SP_basis);

//    sO.system->add_potential(sO.onebody_pot);

    Orbitals* Oxygen = new hydrogenicOrbitals(generalParams, variationalParams);
    systemObjects.SP_basis = Oxygen;
    Fermions system(generalParams, Oxygen);
    system.add_potential(new AtomCore(generalParams));
    system.add_potential(new Coulomb(generalParams));
    systemObjects.system = &system;

    systemObjects.jastrow = new Pade_Jastrow(generalParams, variationalParams);
//    systemObjects.jastrow = new No_Jastrow();
    systemObjects.sample_method = new Importance(generalParams);

    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;

    //Creating the solver objects
    bool silent = true;
    VMC vmc(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w, silent);
    vmc.set_error_estimator(new Blocking(vmcParams.n_c, parParams));

    arma::wall_clock a;
    a.tic();
//    minimizer.minimize();
    vmc.run_method();
    vmc.dump_subsamples();

    if (parParams.is_master) std::cout << "Time: " <<setprecision(3) << fixed << a.toc()/60 << std::endl;


    if (parParams.is_master) cout << "~.* QMC fin *.~" << endl;

#ifdef MPI_ON
    MPI_Finalize();
#endif

    return 0;
}


#endif //QMC2_MAIN
