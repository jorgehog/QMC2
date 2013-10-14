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
#include "Orbitals/NBodyTransform/nbodytransform.h"
//#include "Orbitals/ExpandedBasis/ExpandedBasis.h"

#include "Potential/AtomCore/AtomCore.h"
#include "Potential/Coulomb/Coulomb.h"
#include "Potential/DiAtomCore/DiAtomCore.h"
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


    minimizerParams.SGDsamples = 3000;

    dmcParams.therm=1000;
    dmcParams.n_c = 2000;
    vmcParams.n_c = 1E6;
    dmcParams.dt = 0.00005;

    //    variationalParams.alpha = 1.35618;
    //    variationalParams.beta =  0.28;
    //    generalParams.R = 1.4;

    variationalParams.alpha = 1.1451;
    variationalParams.beta =  0.31;

    //Setting up parallel parameters
    initMPI(parParams, argc, argv);

    dmcParams.n_w = 200*parParams.n_nodes;
    minimizerParams.n_c_SGD=100*parParams.n_nodes;

    scaleWithProcs(parParams, generalParams, minimizerParams, vmcParams, dmcParams);

    //Setting up the solver parameters
    generalParams.dim = 3;
    generalParams.n_p = 30;

    double R = 3.0;       //bohr radii
    double theta = 140.0; //deg

    theta *= arma::datum::pi/180;//dmcE: -433.274597 | E_T: -433.357147 | Nw:   957 | 100.0%


    BodyDef Silicon;
    Silicon.n_p_local = 14;
    Silicon.origin << 0 << 0 << 0;

    BodyDef Oxygen1;
    Oxygen1.n_p_local = 8;
    Oxygen1.origin << R << 0 << 0;

    BodyDef Oxygen2;
    Oxygen2.n_p_local = 8;
    Oxygen2.origin << R*cos(theta) << R*sin(theta) << 0;

    std::vector<BodyDef> bodies;
    bodies.push_back(Silicon);
    bodies.push_back(Oxygen1);
    bodies.push_back(Oxygen2);

    NBodyTransform *Molecule = new NBodyTransform(generalParams, variationalParams, ATOMS, bodies);
    systemObjects.SP_basis = Molecule;
    Fermions system(generalParams, Molecule);
    system.add_potential(new MolecularCoulomb(generalParams, Molecule));
    system.add_potential(new Coulomb(generalParams));
    systemObjects.system = &system;

    systemObjects.jastrow = new Pade_Jastrow(generalParams, variationalParams);
    systemObjects.sample_method = new Importance(generalParams);

    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;

    //Creating the solver objects
    bool silent = true;
    VMC vmc(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w, silent);
    ASGD minimizer(&vmc, minimizerParams, parParams, generalParams.runpath);
    DMC dmc(generalParams, dmcParams, systemObjects, parParams, silent);

    vmc.set_error_estimator(new Blocking(vmcParams.n_c, parParams));
    dmc.set_error_estimator(new SimpleVar(parParams));

    //    forceLoop(minimizer, vmc, dmc, &generalParams.R, generalParams.n_p, parParams.is_master);

    arma::wall_clock a;
    a.tic();
//    minimizer.minimize();
    vmc.run_method();
    dmc.initFromVMC(&vmc);
    dmc.run_method();
    if (parParams.is_master) std::cout << "Time: " <<setprecision(3) << fixed << a.toc()/60 << std::endl;

    if (parParams.is_master) cout << "~.* QMC fin *.~" << endl;

#ifdef MPI_ON
    MPI_Finalize();
#endif

    return 0;
}


void forceLoop(ASGD &asgd, VMC &vmc, DMC &dmc, double *R, int n_p, bool is_master){

    using namespace arma;

    double f_dmc, f_vmc, f0;

    int nPoints = 200;
    double rMin = 0.1;
    double rMax = 6;

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

//    dmc.initFromVMC(&vmc);
//    dmc.run_method();

//    f_dmc = forceSampler.extract_mean_of_means();

//    f0 = 2*f_dmc - f_vmc;
//    fVec(0) = f0;

    fVec(0) = f_vmc;

    fVec.save("/tmp/binaryArmaVec.arma");

    for (int i = 1; i < nPoints; ++i) {

        *R = rVec(i);

        if (is_master) std::cout << "\nLooping R = " << *R << "\n" << std::endl;

        asgd.minimize(false);
        vmc.run_method(false);

        f_vmc = forceSampler.extract_mean();

//        dmc.run_method(false);

//        f_dmc = forceSampler.extract_mean_of_means();

//        f0 = 2*f_dmc - f_vmc;
//        fVec(i) = f0;

        fVec(i) = f_vmc;

        fVec.save("/tmp/binaryArmaVec.arma");

    }
}

#endif //QMC2_MAIN
