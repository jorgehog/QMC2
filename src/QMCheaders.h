/* 
 * File:   QMCheaders.h
 * Author: jorgehog
 *
 * Created on 26. juni 2012, 15:19
 */

#ifndef QMCHEADERS_H
#define	QMCHEADERS_H

//#define ARMA_NO_DEBUG

#include <armadillo>
#include <stdlib.h>
#include <iostream> 
#include <fstream>
#include <math.h>
#include <vector>

struct DMCparams {
    int n_c, therm, n_w, n_b;
    double dt, E_T;
    bool dist_in;
    std::string dist_in_path;

};

struct VMCparams {
    int n_c;
    double dt;

};

struct VariationalParams {
    double alpha, beta;

};

struct GeneralParams {
    int n_p, dim;
    long random_seed;
    double D;

    double h, w;

    int numprocs, myrank;

    bool parallell;

    bool doMIN;
    bool doVMC;
    bool doDMC;

    bool use_jastrow;
    bool use_coulomb;

    std::string system;
    std::string sampling;

};

struct MinimizerParams {
    double max_step;
    double f_max;
    double f_min;
    double omega;
    double A;
    double a;
    int SGDsamples;
    int n_walkers;
    int thermalization;
    int n_cm;
    int n_c_SGD;
    arma::rowvec alpha;
    arma::rowvec beta;

};

struct OutputParams {
    bool dist_out;
    bool blocking_out;
    bool dmc_out;
    std::string outputSuffix;
    std::string outputPath;

};

class Orbitals;
class Potential;
class System;
class Sampling;
class Jastrow;

struct SystemObjects {
    Orbitals* SP_basis;
    Potential* onebody_pot;
    System* SYSTEM;
    Sampling* sample_method;
    Jastrow* jastrow;

};

#include "Walker/Walker.h"

class Minimizer;
class OutputHandler;

#include "BasisFunctions/BasisFunctions.h"
#include "BasisFunctions/semiAlphaHO/semiAlphaHO.h"
#include "BasisFunctions/alphaHO/alphaHO.h"

#include "Orbitals/Orbitals.h"
#include "Orbitals/AlphaHarmonicOscillatorOld/AlphaHarmonicOscillatorOld.h"
#include "Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h"
#include "Orbitals/ExpandedBasis/ExpandedBasis.h"

#include "Potential/Potential.h"
#include "Potential/Coulomb/Coulomb.h"
#include "Potential/Harmonic_osc/Harmonic_osc.h"

#include "Jastrow/Jastrow.h"
#include "Jastrow/Pade_Jastrow/Pade_Jastrow.h"
#include "Jastrow/No_Jastrow/No_Jastrow.h"

#include "System/System.h"
#include "System/Fermions/Fermions.h"

#include "QMC/QMC.h"
#include "QMC/VMC/VMC.h"
#include "QMC/DMC/DMC.h"

#include "Diffusion/Diffusion.h"
#include "Diffusion/Simple/Simple.h"
#include "Diffusion/Fokker_Planck/Fokker_Planck.h"

#include "Sampling/Sampling.h"
#include "Sampling/Brute_Force/Brute_Force.h"
#include "Sampling/Importance/Importance.h"

#include "OutputHandler/OutputHandler.h"
#include "OutputHandler/Distribution/Distribution.h"
#include "OutputHandler/BlockingData/BlockingData.h"
#include "OutputHandler/stdoutDMC/stdoutDMC.h"

#include "Minimizer/Minimizer.h"
#include "Minimizer/ASGD/ASGD.h"


#endif	/* QMCHEADERS_H */

