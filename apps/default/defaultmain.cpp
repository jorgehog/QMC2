#include <QMC2.h>

/*
 *
 */

using namespace QMC2;

namespace QMC2 {

namespace SYSTEMS {

enum {

    QDots,
    DoubleWell,
    Atoms,
    Diatom,
    QDots3D

};
}

namespace SAMPLING {

enum {
    IS,
    BF
};
}

}

struct MainFileParams {

    bool doMIN = false;
    bool doVMC = true;
    bool doDMC = false;
    bool do_blocking = false;
    bool use_jastrow = true;
    bool use_coulomb = true;

    int system = QMC2::SYSTEMS::QDots;
    int sampling = QMC2::SAMPLING::IS;

};



//! Function for parsing the command line for parameters.
/*!
 * The command line parameters are set up by the Python script runQMC.
 * Compiling with a different main file is designed to be easy.
 */
void parseCML(int argc, char** argv,
              MainFileParams & mP,
              VMCparams & vmcParams,
              DMCparams & dmcParams,
              VariationalParams & variationalParams,
              GeneralParams & generalParams,
              MinimizerParams & minimizerParams,
              ParParams & parParams);


//! Function for setting up the system objects.
void selectSystem(MainFileParams &mP, GeneralParams & gP,
                  SystemObjects & sO,
                  VariationalParams & vP,
                  ParParams & pp);

void blocking(int argc, char**argv, struct ParParams & parParams);

void dist(int argc, char**argv, struct ParParams & parParams);

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

    //Initializing solver objects
    Minimizer* minimizer;
    VMC* vmc;
    DMC* dmc;

    //Setting up parallel parameters
    initMPI(parParams, argc, argv);

    //THIS IS CONTROLLED BY PYTHON SCRIPT
    blocking(argc, argv, parParams);
    dist(argc, argv, parParams);

    struct MainFileParams mainParams;

    parseCML(argc, argv,
             mainParams,
             vmcParams,
             dmcParams,
             variationalParams,
             generalParams,
             minimizerParams,
             parParams);
    //END OF PYTHON CONTROLLED ROUTINES

    arma::wall_clock t;
    scaleWithProcs(parParams, generalParams, minimizerParams, vmcParams, dmcParams);

    //Deprecate this
    selectSystem(mainParams, generalParams, systemObjects, variationalParams, parParams);


    if (mainParams.doVMC) {

        vmc = new VMC(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w);

        if (mainParams.doMIN) {

            minimizer = new ASGD(vmc, minimizerParams, parParams, generalParams.runpath);

            if (parParams.is_master) t.tic();
            minimizer->minimize();
            if (parParams.is_master) cout << "---Minimization time: " << t.toc() << " s---\n" << endl;

        }


        if (mainParams.do_blocking) {

            ErrorEstimator* blocking = new Blocking(vmcParams.n_c,
                                                    parParams,
                                                    "blocking_VMC_out",
                                                    generalParams.runpath);

            vmc->set_error_estimator(blocking);
        } else {
            vmc->set_error_estimator(new SimpleVar(parParams));
        }


        if (parParams.is_master) t.tic();
        vmc->run_method();
        if (parParams.is_master) cout << "---VMC time: " << t.toc() << " s---\n" << endl;

    } else {
        vmc = NULL;
    }

    if (mainParams.doDMC) {

        dmc = new DMC(generalParams, dmcParams, systemObjects, parParams, vmc);

        if (mainParams.do_blocking) {

            ErrorEstimator* blocking = new Blocking(dmcParams.n_c,
                                                    parParams,
                                                    "blocking_DMC_out",
                                                    generalParams.runpath);

            dmc->set_error_estimator(blocking);
        } else {
            dmc->set_error_estimator(new SimpleVar(parParams));
        }

        if (parParams.is_master) t.tic();
        dmc->run_method();
        if (parParams.is_master) cout << "---DMC time: " << t.toc() << " s---\n" << endl;

    }

    if (parParams.is_master) cout << "~.* QMC fin *.~" << endl;

#ifdef MPI_ON
    MPI_Finalize();
#endif

    return 0;
}



void selectSystem(MainFileParams &mP,
                  GeneralParams & gP,
                  SystemObjects & sO,
                  VariationalParams & vP,
                  ParParams & pp) {

    using namespace QMC2;

    System* system;
    Potential* onebody_pot;

    switch (mP.system) {
    case SYSTEMS::QDots:
        gP.dim = 2;

        sO.SP_basis = new AlphaHarmonicOscillator(gP, vP);

        onebody_pot = new Harmonic_osc(gP);

        sO.system = new Fermions(gP, sO.SP_basis);
        sO.system->add_potential(onebody_pot);

        break;

    case SYSTEMS::QDots3D:

        gP.dim = 3;

        sO.SP_basis = new AlphaHarmonicOscillator(gP, vP);

        onebody_pot = new Harmonic_osc(gP);

        sO.system = new Fermions(gP, sO.SP_basis);
        sO.system->add_potential(onebody_pot);

        break;

    case SYSTEMS::Atoms:

        gP.dim = 3;

        sO.SP_basis = new hydrogenicOrbitals(gP, vP);

        onebody_pot = new AtomCore(gP);

        sO.system = new Fermions(gP, sO.SP_basis);

        sO.system->add_potential(onebody_pot);

        break;

    case SYSTEMS::Diatom:

        gP.dim = 3;

//        sO.SP_basis = new DiTransform(gP, vP, ATOMS);

        system = new Fermions(gP, sO.SP_basis);
//        system->add_potential(new DiAtomCore(gP));

        sO.system = system;

        break;

    case SYSTEMS::DoubleWell:

        gP.dim = 2;

//        sO.SP_basis = new DiTransform(gP, vP, QDOTS);

        system = new Fermions(gP, sO.SP_basis);
        system->add_potential(new DoubleWell(gP, 1));

        sO.system = system;

        break;

    default:

        if (pp.is_master) std::cout << "unknown system " << system << std::endl;
        if (pp.is_master) exit(1);

        break;
    }


    switch (mP.sampling) {
    case SAMPLING::IS:

        sO.sample_method = new Importance(gP);

        break;

    case SAMPLING::BF:

        sO.sample_method = new Brute_Force(gP);

        break;

    default:

        if (pp.is_master) std::cout << "unknown sampling method" << std::endl;
        if (pp.is_master) exit(1);

        break;
    }



    if (mP.use_jastrow) {
        sO.jastrow = new Pade_Jastrow(gP, vP);

    } else {
        sO.jastrow = new No_Jastrow();
    }



    if (mP.use_coulomb) {
        sO.system->add_potential(new Coulomb(gP));
    }


}


void parseCML(int argc, char** argv,
              MainFileParams & mP,
              VMCparams & vmcParams,
              DMCparams & dmcParams,
              VariationalParams & variationalParams,
              GeneralParams & generalParams,
              MinimizerParams & minimizerParams,
              ParParams & parParams) {


    int n_args = 32;

    //Seting values if not flagged default (controlled by Python)
    std::string def = "def";

    if (argc == n_args) {


        if (def.compare(argv[1]) != 0) generalParams.runpath = argv[1];
        if (def.compare(argv[2]) != 0) generalParams.n_p = atoi(argv[2]);
        if (def.compare(argv[3]) != 0) generalParams.dim = atoi(argv[3]);
        if (def.compare(argv[4]) != 0) generalParams.systemConstant = atof(argv[4]);
//        if (def.compare(argv[5]) != 0) generalParams.R = atof(argv[5]);



        if (def.compare(argv[6]) != 0) generalParams.random_seed = atoi(argv[6]);



        if (def.compare(argv[7]) != 0) mP.doMIN = (bool)atoi(argv[7]);
        if (def.compare(argv[8]) != 0) mP.doVMC = (bool)atoi(argv[8]);
        if (def.compare(argv[9]) != 0) mP.doDMC = (bool)atoi(argv[9]);




        if (def.compare(argv[10]) != 0) mP.use_coulomb = (bool)atoi(argv[10]);
        if (def.compare(argv[11]) != 0) mP.use_jastrow = (bool)atoi(argv[11]);
        if (def.compare(argv[12]) != 0) mP.do_blocking = (bool)atoi(argv[12]);


        if (def.compare(argv[13]) != 0) mP.sampling = atoi(argv[13]);
        if (def.compare(argv[14]) != 0) mP.system = atoi(argv[14]);


        if (def.compare(argv[15]) != 0) generalParams.deadlock_x = atof(argv[15]);




        if (def.compare(argv[16]) != 0) vmcParams.n_c = atoi(argv[16]);
        if (def.compare(argv[17]) != 0) vmcParams.dt = atof(argv[17]);




        if (def.compare(argv[18]) != 0) dmcParams.dt = atof(argv[18]);
        if (def.compare(argv[19]) != 0) dmcParams.n_b = atoi(argv[19]);
        if (def.compare(argv[20]) != 0) dmcParams.n_w = atoi(argv[20]);
        if (def.compare(argv[21]) != 0) dmcParams.n_c = atoi(argv[21]);
        if (def.compare(argv[22]) != 0) dmcParams.therm = atoi(argv[22]);




        if (def.compare(argv[23]) != 0) minimizerParams.SGDsamples = atoi(argv[23]);
        if (def.compare(argv[24]) != 0) minimizerParams.n_w = atoi(argv[24]);
        if (def.compare(argv[25]) != 0) minimizerParams.therm = atoi(argv[25]);
        if (def.compare(argv[26]) != 0) minimizerParams.n_c_SGD = atoi(argv[26]);
        if (def.compare(argv[27]) != 0) minimizerParams.max_step = atof(argv[27]);
        if (def.compare(argv[28]) != 0) minimizerParams.alpha = arma::zeros(1, 1) + atof(argv[28]);
        if (def.compare(argv[29]) != 0) minimizerParams.beta = arma::zeros(1, 1) + atof(argv[29]);




        if (def.compare(argv[30]) != 0) variationalParams.alpha = atof(argv[30]);
        if (def.compare(argv[31]) != 0) variationalParams.beta = atof(argv[31]);

        int vmc_dt_loc = 17;
        int deadlock_loc = 15;


        if (def.compare(argv[vmc_dt_loc]) == 0) {
            if (mP.sampling == QMC2::SAMPLING::IS) {
                vmcParams.dt = 0.005;
            } else {
                vmcParams.dt = 0.5;
            }
        }

        if (def.compare(argv[deadlock_loc]) != 0)
            generalParams.deadlock = true;

    } else {
        std::cout << "insufficient CML arguments. Attempting execution..." << std::endl;
    }


    if (!(mP.doVMC) && mP.doMIN) {
        //we initialize a VMC run at the end either way. Cheap verification.
        mP.doVMC = 1;
    }

    //Check consitency in chosen cycle numbers etc. given node number.
    if (parParams.is_master) {
        if (parParams.parallel) {
            if (mP.doMIN) {
                if (minimizerParams.n_c_SGD >= parParams.n_nodes) {

                } else {
                    std::cout << "n_c_SGD=" << parParams.n_nodes << " is too low for n_nodes=" << parParams.n_nodes << std::endl;
                    std::cout << "aborting.";
                    exit(1);
                }
            }


        }

        if (mP.doDMC) {
            if (dmcParams.n_w < parParams.n_nodes) {
                std::cout << "n_w = " << dmcParams.n_w << " is insufficient for ";
                std::cout << parParams.n_nodes << " nodes." << std::endl;
            }


            if (mP.doVMC) {

                if (dmcParams.n_w > vmcParams.n_c) {
                    std::cout << "Unsufficient VMC cycles store all walkers." << std::endl;
                    std::cout << "For n_w=" << dmcParams.n_w << "the minimum is n_c=" << dmcParams.n_w << std::endl;
                    exit(1);
                }

            }
        }

    }



    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;

    bool initOut = false;
    if (initOut) {
        if (parParams.is_master) {
            std::cout << n_args << " =? " << argc << std::endl;
            std::cout << " 4" << " generalParams.runpath           " << " = " << generalParams.runpath << "   " << argv[4] << std::endl;
            std::cout << " 5" << " generalParams.n_p               " << " = " << generalParams.n_p << "   " << argv[5] << std::endl;
            std::cout << " 6" << " generalParams.dim               " << " = " << generalParams.dim << "   " << argv[6] << std::endl;
            std::cout << " 7" << " generalParams.systemConstant    " << " = " << generalParams.systemConstant << "   " << argv[7] << std::endl;
//            std::cout << " 8" << " generalParams.R                 " << " = " << generalParams.R << "   " << argv[8] << std::endl;
            std::cout << " 9" << " generalParams.random_seed       " << " = " << generalParams.random_seed << "   " << argv[9] << std::endl;
            std::cout << "10" << " generalParams.doMIN             " << " = " << mP.doMIN << "   " << argv[10] << std::endl;
            std::cout << "11" << " generalParams.doVMC             " << " = " << mP.doVMC << "   " << argv[11] << std::endl;
            std::cout << "12" << " generalParams.doDMC             " << " = " << mP.doDMC << "   " << argv[12] << std::endl;
            std::cout << "13" << " generalParams.use_coulomb       " << " = " << mP.use_coulomb << "   " << argv[13] << std::endl;
            std::cout << "14" << " generalParams.use_jastrow       " << " = " << mP.use_jastrow << "   " << argv[14] << std::endl;
            std::cout << "15" << " generalParams.do_blocking       " << " = " << mP.do_blocking << "   " << argv[15] << std::endl;
            std::cout << "16" << " generalParams.sampling          " << " = " << mP.sampling << "   " << argv[16] << std::endl;
            std::cout << "17" << " generalParams.system            " << " = " << mP.system << "   " << argv[17] << std::endl;
            std::cout << "18" << " vmcParams.n_c                   " << " = " << vmcParams.n_c << "   " << argv[18] << std::endl;
            std::cout << "19" << " vmcParams.dt                    " << " = " << vmcParams.dt << "   " << argv[19] << std::endl;
            std::cout << "20" << " dmcParams.dt                    " << " = " << dmcParams.dt << "   " << argv[20] << std::endl;
            std::cout << "21" << " dmcParams.n_b                   " << " = " << dmcParams.n_b << "   " << argv[21] << std::endl;
            std::cout << "22" << " dmcParams.n_w                   " << " = " << dmcParams.n_w << "   " << argv[22] << std::endl;
            std::cout << "23" << " dmcParams.n_c                   " << " = " << dmcParams.n_c << "   " << argv[23] << std::endl;
            std::cout << "24" << " dmcParams.therm                 " << " = " << dmcParams.therm << "   " << argv[24] << std::endl;
            std::cout << "25" << " minimizerParams.SGDsamples      " << " = " << minimizerParams.SGDsamples << "   " << argv[25] << std::endl;
            std::cout << "26" << " minimizerParams.n_w             " << " = " << minimizerParams.n_w << "   " << argv[26] << std::endl;
            std::cout << "27" << " minimizerParams.therm           " << " = " << minimizerParams.therm << "   " << argv[27] << std::endl;
            std::cout << "28" << " minimizerParams.n_c_SGD         " << " = " << minimizerParams.n_c_SGD << "   " << argv[28] << std::endl;
            std::cout << "29" << " minimizerParams.max_step        " << " = " << minimizerParams.max_step << "   " << argv[29] << std::endl;
            std::cout << "30" << " minimizerParams.alpha           " << " = " << minimizerParams.alpha << "   " << argv[30] << std::endl;
            std::cout << "31" << " minimizerParams.beta            " << " = " << minimizerParams.beta << "   " << argv[31] << std::endl;
            std::cout << "32" << " variationalParams.alpha         " << " = " << variationalParams.alpha << "   " << argv[32] << std::endl;
            std::cout << "33" << " variationalParams.beta          " << " = " << variationalParams.beta << "   " << argv[33] << std::endl;
        }
#ifdef MPI_ON
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
#endif

        exit(0);
    }



#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//Test for rerunning blocking
//argv = [x-name, "reblock", filename, path, #blocks, max_block_size, min_block_size]?
void blocking(int argc, char** argv, struct ParParams & parParams){

    using namespace std;

    string rerun_blocking = "reblock";
    if ((argc == 7) && (rerun_blocking.compare(argv[1]) == 0)) {
        ErrorEstimator* reblock = new Blocking(0, parParams,
                                               (string) argv[2],
                (string) argv[3],
                atoi(argv[4]),
                atoi(argv[5]),
                atoi(argv[6]),
                true);
        double error = reblock->estimate_error();
        if (parParams.is_master) cout << "Estimated Error: " << error << endl;
        if (parParams.is_master) cout << "Finished Error Recalculation" << endl;

        reblock->finalize();

#ifdef MPI_ON
        MPI_Finalize();
#endif
        exit(0);
    }
}
//

//Test for rerunning distribution
//argv = [x-name, "redist", n_p, path, name, N, bin_edge?
void dist(int argc, char** argv, struct ParParams & parParams){

    using namespace std;

    string rerun_dist = "redist";
    if ((argc == 7) && (rerun_dist.compare(argv[1]) == 0)) {
        Distribution* redist = new Distribution(parParams, argv[3], argv[4]);

        redist->rerun(atoi(argv[2]), atoi(argv[5]), atof(argv[6]));

#ifdef MPI_ON
        MPI_Finalize();
#endif
        exit(0);
    }
}
