/* 
 * File:   QMCmain.cpp
 * Author: jorgehog
 *
 * Created on 13. april 2012, 17:04
 */

#include "QMCheaders.h"
//#include "HartreeFock/HartreeFock.h"

/*
 * 
 */

//! Function for parsing the command line for parameters.
/*!
 * The command line parameters are set up by the Python script runQMC.
 * Compiling with a different main file is designed to be easy.
 */
void parseCML(int, char**,
        VMCparams &,
        DMCparams &,
        VariationalParams &,
        GeneralParams &,
        MinimizerParams &,
        OutputParams &,
        ParParams &
        );

//! Function for setting up the system objects.
void selectSystem(GeneralParams & gP,
        SystemObjects & sO,
        VariationalParams & vP,
        ParParams & pp);

int main(int argc, char** argv) {

    using namespace std;

    //Setting up parallel parameters
    struct ParParams parParams;

#ifdef MPI_ON
    int node, n_nodes;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);

    parParams.n_nodes = n_nodes;
    parParams.node = node;
    parParams.parallel = (parParams.n_nodes > 1);
    parParams.is_master = (parParams.node == 0);
#else
    parParams.parallel = false;
    parParams.node = 0;
    parParams.n_nodes = 1;
    parParams.is_master = true;
#endif

    arma::wall_clock t;

    //Initialize structs for holding constructor arguments.
    struct VMCparams vmcParams;
    struct DMCparams dmcParams;
    struct VariationalParams variationalParams;
    struct GeneralParams generalParams;
    struct MinimizerParams minimizerParams;
    struct OutputParams outputParams;
    struct SystemObjects systemObjects;


    //Test for rerunning blocking
    //argv = [x-name, "reblock", filename, path, #blocks, max_block_size, min_block_size]?
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
        return 0;
    }
    //

    //Test for rerunning distribution
    //argv = [x-name, "redist", n_p, path, name, N, bin_edge?
    string rerun_dist = "redist";
    if (rerun_dist.compare(argv[1]) == 0) {
        Distribution* redist = new Distribution(parParams, argv[3], argv[4]);
        
        if (argc == 7) {
            redist->rerun(atoi(argv[2]), atoi(argv[5]), atof(argv[6]));
        } 
        
#ifdef MPI_ON
        MPI_Finalize();
#endif
        return 0;
    }
    
//    generalParams.n_p=2;
//    generalParams.dim=3;
//    generalParams.systemConstant = generalParams.n_p;
//    variationalParams.alpha = 0.6425;
//    variationalParams.beta = 0.28;
//    double* R;
//    *R = 1.4;
//    
//    Orbitals* 
    //

    //    generalParams.n_p = 2;
    //    generalParams.dim = 3;
    //    generalParams.systemConstant = generalParams.n_p;
    //    variationalParams.alpha = 1;
    //    srand(time(NULL));
    //    Orbitals* sp = new hydrogenicOrbitals(generalParams, variationalParams);
    //    HartreeFock* hf = new HartreeFock(4, sp);
    //    hf->run_method();
    //    
    //    exit(0);

    parseCML(argc, argv,
            vmcParams,
            dmcParams,
            variationalParams,
            generalParams,
            minimizerParams,
            outputParams,
            parParams);

    selectSystem(generalParams, systemObjects, variationalParams, parParams);

    Minimizer* minimizer;
    VMC* vmc;
    DMC* dmc;

    std::stringstream name;

    if (generalParams.doVMC) {

        vmc = new VMC(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w, outputParams.dist_out);

        if (generalParams.doMIN) {

            minimizer = new ASGD(vmc, minimizerParams, parParams);

            if (outputParams.ASGD_out && parParams.is_master) {
                OutputHandler* ASGDout = new stdoutASGD(generalParams.runpath);
                minimizer->add_output(ASGDout);
            }


            int N = minimizerParams.alpha.n_elem + minimizerParams.beta.n_rows * generalParams.use_jastrow;
            for (int i = 0; i < N; i++) {
                if (generalParams.do_blocking) {

                    ErrorEstimator* blocking = new Blocking(minimizerParams.SGDsamples,
                            parParams,
                            "blocking_MIN_out" + TOSTR(i),
                            generalParams.runpath);

                    minimizer->add_error_estimator(blocking);
                } else {
                    minimizer->add_error_estimator(new SimpleVar(parParams));
                }
            }

            if (parParams.is_master) t.tic();
            minimizer->minimize();
            if (parParams.is_master) cout << "---Minimization time: " << t.toc() << " s---\n" << endl;
        }



        if (outputParams.dist_out) {

            name << generalParams.system << generalParams.n_p << "c" << generalParams.systemConstant;
            name << "vmc";

            OutputHandler* dist_vmc = new Distribution(parParams, generalParams.runpath, name.str());
            vmc->add_output(dist_vmc);

            name.str(std::string());
            name.clear();

        }


        if (generalParams.do_blocking) {

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

    }

    if (generalParams.doDMC) {

        dmc = new DMC(generalParams, dmcParams, systemObjects, parParams, vmc, outputParams.dist_out);

        if (outputParams.dmc_out && parParams.is_master) {

            OutputHandler* DMCout = new stdoutDMC(generalParams.runpath);
            dmc->add_output(DMCout);
        }

        if (outputParams.dist_out) {

            name << generalParams.system << generalParams.n_p << "c" << generalParams.systemConstant;
            name << "dmc";

            OutputHandler* dist_dmc = new Distribution(parParams, generalParams.runpath, name.str());
            dmc->add_output(dist_dmc);

            name.str(std::string());
            name.clear();

        }

        //        int DMCerrorN = dmcParams.n_c;
        if (generalParams.do_blocking) {

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

void parseCML(int argc, char** argv,
        VMCparams & vmcParams,
        DMCparams & dmcParams,
        VariationalParams & variationalParams,
        GeneralParams & generalParams,
        MinimizerParams & minimizerParams,
        OutputParams & outputParams,
        ParParams & parParams) {


    int n_args = 33;

    //Default values:

    generalParams.n_p = 2;
    generalParams.dim = 2;
    generalParams.systemConstant = 1;
    generalParams.random_seed = (seed_type) time(NULL);
    generalParams.doMIN = false;
    generalParams.doVMC = false;
    generalParams.doDMC = false;

    generalParams.use_coulomb = true;
    generalParams.use_jastrow = true;

    generalParams.sampling = "IS";
    generalParams.do_blocking = false;
    generalParams.system = "QDots";


    outputParams.dist_out = true;
    outputParams.dmc_out = true;
    outputParams.ASGD_out = true;
    generalParams.runpath = "/home/jorgmeister/scratch/QMC_SCRATCH/";

    vmcParams.n_c = 1E6;
    if (generalParams.sampling == "IS") {
        vmcParams.dt = 0.005;
    } else {
        vmcParams.dt = 0.5;
    }


    dmcParams.dt = 0.001;
    dmcParams.n_b = 100;
    dmcParams.n_c = 1000;
    dmcParams.n_w = 1000;
    dmcParams.therm = 1000;

    if (argc == 1) {
        variationalParams.alpha = 0.987;
        variationalParams.beta = 0.398;
    }


    //defauls ASGD parameters
    minimizerParams.max_step = 0.1;
    minimizerParams.f_max = 1.0;
    minimizerParams.f_min = -0.5;
    minimizerParams.omega = 0.8;
    minimizerParams.A = 60;
    minimizerParams.a = 0.5;
    minimizerParams.n_w = 10;
    minimizerParams.therm = 10000;
    minimizerParams.n_c = 1000;

    minimizerParams.SGDsamples = 2000;
    minimizerParams.n_c_SGD = 400;
    minimizerParams.alpha = arma::zeros(1, 1) + 0.5;
    minimizerParams.beta = arma::zeros(1, 1) + 0.5;


    //Seting values if not flagged default (controlled by Python)
    std::string def = "def";

    if (argc == n_args) {
        if (def.compare(argv[1]) != 0) outputParams.dist_out = (bool)atoi(argv[1]);
        if (def.compare(argv[2]) != 0) outputParams.dmc_out = (bool)atoi(argv[2]);
        if (def.compare(argv[3]) != 0) outputParams.ASGD_out = (bool)atoi(argv[3]);



        if (def.compare(argv[4]) != 0) generalParams.runpath = argv[4];
        if (def.compare(argv[5]) != 0) generalParams.n_p = atoi(argv[5]);
        if (def.compare(argv[6]) != 0) generalParams.dim = atoi(argv[6]);
        if (def.compare(argv[7]) != 0) generalParams.systemConstant = atof(argv[7]);


        if (def.compare(argv[8]) != 0) generalParams.random_seed = atoi(argv[8]);


        if (def.compare(argv[9]) != 0) generalParams.doMIN = (bool)atoi(argv[9]);
        if (def.compare(argv[10]) != 0) generalParams.doVMC = (bool)atoi(argv[10]);
        if (def.compare(argv[11]) != 0) generalParams.doDMC = (bool)atoi(argv[11]);



        if (def.compare(argv[12]) != 0) generalParams.use_coulomb = (bool)atoi(argv[12]);
        if (def.compare(argv[13]) != 0) generalParams.use_jastrow = (bool)atoi(argv[13]);
        if (def.compare(argv[14]) != 0) generalParams.do_blocking = (bool)atoi(argv[14]);



        if (def.compare(argv[15]) != 0) generalParams.sampling = argv[15];
        if (def.compare(argv[16]) != 0) generalParams.system = argv[16];



        if (def.compare(argv[17]) != 0) vmcParams.n_c = atoi(argv[17]);
        if (def.compare(argv[18]) != 0) vmcParams.dt = atof(argv[18]);



        if (def.compare(argv[19]) != 0) dmcParams.dt = atof(argv[19]);
        if (def.compare(argv[20]) != 0) dmcParams.n_b = atoi(argv[20]);
        if (def.compare(argv[21]) != 0) dmcParams.n_w = atoi(argv[21]);
        if (def.compare(argv[22]) != 0) dmcParams.n_c = atoi(argv[22]);
        if (def.compare(argv[23]) != 0) dmcParams.therm = atoi(argv[23]);



        if (def.compare(argv[24]) != 0) minimizerParams.SGDsamples = atoi(argv[24]);
        if (def.compare(argv[25]) != 0) minimizerParams.n_w = atoi(argv[25]);
        if (def.compare(argv[26]) != 0) minimizerParams.therm = atoi(argv[26]);
        if (def.compare(argv[27]) != 0) minimizerParams.n_c_SGD = atoi(argv[27]);
        if (def.compare(argv[28]) != 0) minimizerParams.max_step = atof(argv[28]);
        if (def.compare(argv[29]) != 0) minimizerParams.alpha = arma::zeros(1, 1) + atof(argv[29]);
        if (def.compare(argv[30]) != 0) minimizerParams.beta = arma::zeros(1, 1) + atof(argv[30]);



        if (def.compare(argv[31]) != 0) variationalParams.alpha = atof(argv[31]);
        if (def.compare(argv[32]) != 0) variationalParams.beta = atof(argv[32]);


        int vmc_dt_loc = 18;

        if (def.compare(argv[vmc_dt_loc]) == 0) {
            if (generalParams.sampling == "IS") {
                vmcParams.dt = 0.005;
            } else {
                vmcParams.dt = 0.5;
            }
        }

    } else {
        std::cout << "insufficient CML arguments. Attempting execution..." << std::endl;
    }


    if (!(generalParams.doVMC) && generalParams.doMIN) {
        //we initialize a VMC run at the end either way. Cheap verification.
        generalParams.doVMC = 1;
    }

    //Check consitency in chosen cycle numbers etc. given node number.
    if (parParams.is_master) {
        if (parParams.parallel) {
            if (generalParams.doMIN) {
                if (minimizerParams.n_c_SGD >= parParams.n_nodes) {

                } else {
                    std::cout << "n_c_SGD=" << parParams.n_nodes << " is too low for n_nodes=" << parParams.n_nodes << std::endl;
                    std::cout << "aborting.";
                    exit(1);
                }
            }


        }

        if (generalParams.doDMC) {
            if (dmcParams.n_w < parParams.n_nodes) {
                std::cout << "n_w = " << dmcParams.n_w << " is insufficient for ";
                std::cout << parParams.n_nodes << " nodes." << std::endl;
            }


            if (generalParams.doVMC) {

                if (dmcParams.n_w > vmcParams.n_c) {
                    std::cout << "Unsufficient VMC cycles store all walkers." << std::endl;
                    std::cout << "For n_w=" << dmcParams.n_w << "the minimum is n_c=" << dmcParams.n_w << std::endl;
                    exit(1);
                }

            }
        }

    }

    generalParams.random_seed -= parParams.node;
    minimizerParams.n_c_SGD /= parParams.n_nodes;

    vmcParams.n_c /= parParams.n_nodes;
    dmcParams.n_w /= parParams.n_nodes;

    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;

    bool initOut = false;
    if (initOut) {
        if (parParams.is_master) {
            if (n_args == argc) {

                std::cout << n_args << " =? " << argc << std::endl;
                std::cout << " 1" << " outputParams.dist_out           " << " = " << outputParams.dist_out << "   " << argv[1] << std::endl;
                std::cout << " 2" << " outputParams.dmc_out            " << " = " << outputParams.dmc_out << "   " << argv[2] << std::endl;
                std::cout << " 3" << " outputParams.ASGD_out           " << " = " << outputParams.ASGD_out << "   " << argv[3] << std::endl;
                std::cout << " 4" << " generalParams.runpath           " << " = " << generalParams.runpath << "   " << argv[4] << std::endl;
                std::cout << " 5" << " generalParams.n_p               " << " = " << generalParams.n_p << "   " << argv[5] << std::endl;
                std::cout << " 6" << " generalParams.dim               " << " = " << generalParams.dim << "   " << argv[6] << std::endl;
                std::cout << " 7" << " generalParams.systemConstant    " << " = " << generalParams.systemConstant << "   " << argv[7] << std::endl;
                std::cout << " 8" << " generalParams.random_seed       " << " = " << generalParams.random_seed << "   " << argv[8] << std::endl;
                std::cout << " 9" << " generalParams.doMIN             " << " = " << generalParams.doMIN << "   " << argv[9] << std::endl;
                std::cout << "10" << " generalParams.doVMC             " << " = " << generalParams.doVMC << "   " << argv[10] << std::endl;
                std::cout << "11" << " generalParams.doDMC             " << " = " << generalParams.doDMC << "   " << argv[11] << std::endl;
                std::cout << "12" << " generalParams.use_coulomb       " << " = " << generalParams.use_coulomb << "   " << argv[12] << std::endl;
                std::cout << "13" << " generalParams.use_jastrow       " << " = " << generalParams.use_jastrow << "   " << argv[13] << std::endl;
                std::cout << "14" << " generalParams.do_blocking       " << " = " << generalParams.do_blocking << "   " << argv[14] << std::endl;
                std::cout << "15" << " generalParams.sampling          " << " = " << generalParams.sampling << "   " << argv[15] << std::endl;
                std::cout << "16" << " generalParams.system            " << " = " << generalParams.system << "   " << argv[16] << std::endl;
                std::cout << "17" << " vmcParams.n_c                   " << " = " << vmcParams.n_c << "   " << argv[17] << std::endl;
                std::cout << "18" << " vmcParams.dt                    " << " = " << vmcParams.dt << "   " << argv[18] << std::endl;
                std::cout << "19" << " dmcParams.dt                    " << " = " << dmcParams.dt << "   " << argv[19] << std::endl;
                std::cout << "20" << " dmcParams.n_b                   " << " = " << dmcParams.n_b << "   " << argv[20] << std::endl;
                std::cout << "21" << " dmcParams.n_w                   " << " = " << dmcParams.n_w << "   " << argv[21] << std::endl;
                std::cout << "22" << " dmcParams.n_c                   " << " = " << dmcParams.n_c << "   " << argv[22] << std::endl;
                std::cout << "23" << " dmcParams.therm                 " << " = " << dmcParams.therm << "   " << argv[23] << std::endl;
                std::cout << "24" << " minimizerParams.SGDsamples      " << " = " << minimizerParams.SGDsamples << "   " << argv[24] << std::endl;
                std::cout << "25" << " minimizerParams.n_w             " << " = " << minimizerParams.n_w << "   " << argv[25] << std::endl;
                std::cout << "26" << " minimizerParams.therm           " << " = " << minimizerParams.therm << "   " << argv[26] << std::endl;
                std::cout << "27" << " minimizerParams.n_c_SGD         " << " = " << minimizerParams.n_c_SGD << "   " << argv[27] << std::endl;
                std::cout << "28" << " minimizerParams.max_step        " << " = " << minimizerParams.max_step << "   " << argv[28] << std::endl;
                std::cout << "29" << " minimizerParams.alpha           " << " = " << minimizerParams.alpha << "   " << argv[29] << std::endl;
                std::cout << "30" << " minimizerParams.beta            " << " = " << minimizerParams.beta << "   " << argv[30] << std::endl;
                std::cout << "31" << " variationalParams.alpha         " << " = " << variationalParams.alpha << "   " << argv[31] << std::endl;
                std::cout << "32" << " variationalParams.beta          " << " = " << variationalParams.beta << "   " << argv[32] << std::endl;

            } else {
                std::cout << " 1" << " outputParams.dist_out           " << " = " << outputParams.dist_out << std::endl;
                std::cout << " 2" << " outputParams.dmc_out            " << " = " << outputParams.dmc_out << std::endl;
                std::cout << " 3" << " outputParams.ASGD_out           " << " = " << outputParams.ASGD_out << std::endl;
                std::cout << " 4" << " generalParams.runpath           " << " = " << generalParams.runpath << std::endl;
                std::cout << " 5" << " generalParams.n_p               " << " = " << generalParams.n_p << std::endl;
                std::cout << " 6" << " generalParams.dim               " << " = " << generalParams.dim << std::endl;
                std::cout << " 7" << " generalParams.systemConstant    " << " = " << generalParams.systemConstant << std::endl;
                std::cout << " 8" << " generalParams.random_seed       " << " = " << generalParams.random_seed << std::endl;
                std::cout << " 9" << " generalParams.doMIN             " << " = " << generalParams.doMIN << std::endl;
                std::cout << "10" << " generalParams.doVMC             " << " = " << generalParams.doVMC << std::endl;
                std::cout << "11" << " generalParams.doDMC             " << " = " << generalParams.doDMC << std::endl;
                std::cout << "12" << " generalParams.use_coulomb       " << " = " << generalParams.use_coulomb << std::endl;
                std::cout << "13" << " generalParams.use_jastrow       " << " = " << generalParams.use_jastrow << std::endl;
                std::cout << "14" << " generalParams.do_blocking       " << " = " << generalParams.do_blocking << std::endl;
                std::cout << "15" << " generalParams.sampling          " << " = " << generalParams.sampling << std::endl;
                std::cout << "16" << " generalParams.system            " << " = " << generalParams.system << std::endl;
                std::cout << "17" << " vmcParams.n_c                   " << " = " << vmcParams.n_c << std::endl;
                std::cout << "18" << " vmcParams.dt                    " << " = " << vmcParams.dt << std::endl;
                std::cout << "19" << " dmcParams.dt                    " << " = " << dmcParams.dt << std::endl;
                std::cout << "20" << " dmcParams.n_b                   " << " = " << dmcParams.n_b << std::endl;
                std::cout << "21" << " dmcParams.n_w                   " << " = " << dmcParams.n_w << std::endl;
                std::cout << "22" << " dmcParams.n_c                   " << " = " << dmcParams.n_c << std::endl;
                std::cout << "23" << " dmcParams.therm                 " << " = " << dmcParams.therm << std::endl;
                std::cout << "24" << " minimizerParams.SGDsamples      " << " = " << minimizerParams.SGDsamples << std::endl;
                std::cout << "25" << " minimizerParams.n_w             " << " = " << minimizerParams.n_w << std::endl;
                std::cout << "26" << " minimizerParams.therm           " << " = " << minimizerParams.therm << std::endl;
                std::cout << "27" << " minimizerParams.n_c_SGD         " << " = " << minimizerParams.n_c_SGD << std::endl;
                std::cout << "28" << " minimizerParams.alpha           " << " = " << minimizerParams.alpha << std::endl;
                std::cout << "29" << " minimizerParams.beta            " << " = " << minimizerParams.beta << std::endl;
                std::cout << "30" << " variationalParams.alpha         " << " = " << variationalParams.alpha << std::endl;
                std::cout << "31" << " variationalParams.beta          " << " = " << variationalParams.beta << std::endl;
            }
        }
#ifdef MPI_ON
        MPI_Finalize();
#endif

        exit(0);
    }


#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void selectSystem(GeneralParams & gP,
        SystemObjects & sO,
        VariationalParams & vP,
        ParParams & pp) {


    if (gP.sampling == "IS") {
        sO.sample_method = new Importance(gP);

    } else if (gP.sampling == "BF") {
        sO.sample_method = new Brute_Force(gP);

    } else {

        if (pp.is_master) std::cout << "unknown sampling method" << std::endl;
        if (pp.is_master) exit(1);
    }


    if (gP.use_jastrow) {
        sO.jastrow = new Pade_Jastrow(gP, vP);

    } else {
        sO.jastrow = new No_Jastrow();

    }

    if (gP.system == "QDots") {
        sO.SP_basis = new AlphaHarmonicOscillator(gP, vP);

        sO.onebody_pot = new Harmonic_osc(gP);

        sO.SYSTEM = new Fermions(gP, sO.SP_basis);
        sO.SYSTEM->add_potential(sO.onebody_pot);


    } else if (gP.system == "Atoms") {

        sO.SP_basis = new hydrogenicOrbitals(gP, vP, gP.n_p);

        sO.onebody_pot = new AtomCore(gP);

        sO.SYSTEM = new Fermions(gP, sO.SP_basis);

        sO.SYSTEM->add_potential(sO.onebody_pot);


    } else {
        if (pp.is_master) std::cout << "unknown system" << std::endl;
        if (pp.is_master) exit(1);
    }

    if (gP.use_coulomb) {
        sO.SYSTEM->add_potential(new Coulomb(gP));
    }

}
