/* 
 * File:   QMCmain.cpp
 * Author: jorgehog
 *
 * Created on 13. april 2012, 17:04
 */

#include "QMCheaders.h"

/*
 * 
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

    struct VMCparams vmcParams;
    struct DMCparams dmcParams;
    struct VariationalParams variationalParams;
    struct GeneralParams generalParams;
    struct MinimizerParams minimizerParams;
    struct OutputParams outputParams;
    struct SystemObjects systemObjects;



    //Test for rerunning blocking
    //argv = [x-name, "reblock", filename, path, #blocks, max_block_size, min_block_size]?
    std::string rerun_blocking = "reblock";
    if ((argc == 7) && (rerun_blocking.compare(argv[1]) == 0)) {
        ErrorEstimator* reblock = new Blocking(0, parParams,
                (std::string)argv[2],
                (std::string)argv[3],
                atoi(argv[4]),
                atoi(argv[5]),
                atoi(argv[6]),
                true);
        double error = reblock->estimate_error();
        if (parParams.is_master) std::cout << "Estimated Error: " << error << std::endl;
        if (parParams.is_master) std::cout << "Finished Error Recalculation" << std::endl;

        reblock->finalize();

#ifdef MPI_ON
        MPI_Finalize();
#endif
        return 0;
    }
    //

    parseCML(argc, argv,
            vmcParams,
            dmcParams,
            variationalParams,
            generalParams,
            minimizerParams,
            outputParams,
            parParams);

    if (generalParams.sampling == "IS") {
        systemObjects.sample_method = new Importance(generalParams);
    } else if (generalParams.sampling == "BF") {
        systemObjects.sample_method = new Brute_Force(generalParams);
    } else {
        if (parParams.is_master) cout << "unknown sampling method" << endl;
        if (parParams.is_master) exit(1);
    }


    if (generalParams.use_jastrow) {
        systemObjects.jastrow = new Pade_Jastrow(generalParams, variationalParams);
    } else {
        systemObjects.jastrow = new No_Jastrow();
    }

    if (generalParams.system == "QDots") {
        systemObjects.SP_basis = new AlphaHarmonicOscillator(generalParams, variationalParams);

        systemObjects.onebody_pot = new Harmonic_osc(generalParams);

        systemObjects.SYSTEM = new Fermions(generalParams, systemObjects.SP_basis);
        systemObjects.SYSTEM->add_potential(systemObjects.onebody_pot);


    } else {
        if (parParams.is_master) cout << "unknown system" << endl;
        if (parParams.is_master) exit(1);
    }

    if (generalParams.use_coulomb) {
        systemObjects.SYSTEM->add_potential(new Coulomb(generalParams));
    }



    if (generalParams.doVMC || generalParams.doMIN) {

        VMC* vmc = new VMC(generalParams, vmcParams, systemObjects, parParams);
        systemObjects.sample_method->set_dt(vmcParams.dt);



        if (generalParams.doMIN) {

            Minimizer * minimizer = new ASGD(vmc, minimizerParams, parParams);

            if (outputParams.ASGD_out && parParams.is_master) {
                string ASGDoutname = "ASGD_out" + outputParams.outputSuffix;
                OutputHandler* ASGDout = new stdoutASGD(ASGDoutname, outputParams.outputPath);
                minimizer->add_output(ASGDout);
            }


            int N = minimizerParams.alpha.n_elem + minimizerParams.beta.n_rows * generalParams.use_jastrow;

            for (int i = 0; i < N; i++) {
                if (generalParams.estimate_error) {

                    string minBlockname = (string) "blocking_MIN_out" + outputParams.outputSuffix;
                    minBlockname = minBlockname + boost::lexical_cast<std::string > (i);
                    ErrorEstimator* blocking = new Blocking(minimizerParams.SGDsamples, parParams,
                            minBlockname, outputParams.outputPath);
                    minimizer->add_error_estimator(blocking);
                } else {
                    minimizer->add_error_estimator(new SimpleVar(minimizerParams.SGDsamples, parParams));
                }
            }

            if (parParams.is_master) t.tic();
            vmc = minimizer->minimize();
            if (parParams.is_master) cout << "---Minimization time: " << t.toc() << " s---\n" << endl;
        }

        if (generalParams.doVMC) {

            if (outputParams.dist_out) {
                string distname = "dist_out" + outputParams.outputSuffix;
                OutputHandler* dist = new Distribution(parParams, distname, outputParams.outputPath);
                vmc->add_output(dist);

            }


            if (generalParams.estimate_error) {
                string vmcBlockname = (string) "blocking_VMC_out" + outputParams.outputSuffix;
                ErrorEstimator* blocking = new Blocking(vmcParams.n_c, parParams, vmcBlockname, outputParams.outputPath);
                vmc->set_error_estimator(blocking);
            } else {
                vmc->set_error_estimator(new SimpleVar(vmcParams.n_c, parParams));
            }


            if (parParams.is_master) t.tic();
            vmc->run_method();
            if (parParams.is_master) cout << "---VMC time: " << t.toc() << " s---\n" << endl;

            dmcParams.E_T = vmc->get_energy();
        }


    }

    if (generalParams.doDMC) {


        systemObjects.sample_method->set_dt(dmcParams.dt);

        DMC* dmc = new DMC(generalParams, dmcParams, systemObjects, parParams);

        if (outputParams.dmc_out) {
            string dmcOutname = (string) "DMC_out" + outputParams.outputSuffix;
            OutputHandler* DMCout = new stdoutDMC(dmcOutname, outputParams.outputPath,
                    parParams.parallel, parParams.node, parParams.n_nodes);
            dmc->add_output(DMCout);
        }

        int DMCerrorN = dmcParams.n_c * dmcParams.n_b * dmcParams.n_w * DMC::K;
        if (generalParams.estimate_error) {
            string dmcBlockname = (string) "blocking_DMC_out" + outputParams.outputSuffix;
            ErrorEstimator* blocking = new Blocking(DMCerrorN, parParams, dmcBlockname, outputParams.outputPath);
            dmc->set_error_estimator(blocking);
        } else {
            dmc->set_error_estimator(new SimpleVar(DMCerrorN, parParams));
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


    int n_args = 44;

    //Default values:
    //    outputParams.blocking_out = false;
    outputParams.dist_out = true;
    //NEW
    outputParams.dmc_out = false;
    outputParams.ASGD_out = false;
    outputParams.outputSuffix = "";
    outputParams.outputPath = (std::string)"/home/jorgmeister/scratch/QMC_SCRATCH" + (std::string)"/";

    generalParams.n_p = 2;
    generalParams.dim = 2;
    generalParams.w = 1;
//    generalParams.random_seed = -(long)time(NULL);
    generalParams.random_seed = -1355160055;
    generalParams.h = 0.0001;
    generalParams.D = 0.5;

    generalParams.doMIN = argc == 0;
    generalParams.doVMC = argc == 0;
    generalParams.doDMC = argc == 1;

    generalParams.use_coulomb = true;
    generalParams.use_jastrow = true;

    generalParams.sampling = "IS";
    generalParams.estimate_error = false;
    //    generalParams.kinetics_type = "CF";
    generalParams.system = "QDots";


    vmcParams.n_c = 1E6;

    dmcParams.dt = 0.001;
    dmcParams.E_T = 0;
    dmcParams.n_b = 100;
    dmcParams.n_w = 1000;
    dmcParams.n_c = 1000;
    dmcParams.therm = 1000;
    dmcParams.dist_in = outputParams.dist_out & generalParams.doVMC;
    dmcParams.dist_in_path = (std::string)"/home/jorgmeister/scratch/QMC_SCRATCH" + (std::string)"/";

    if (argc == 1) {
        //        variationalParams.alpha = 0.78;
        //        variationalParams.beta = 0.85;
        variationalParams.alpha = 0.987;
        variationalParams.beta = 0.398;
    }

    minimizerParams.max_step = 0.1;
    minimizerParams.f_max = 1.0;
    minimizerParams.f_min = -0.5;
    minimizerParams.omega = 0.8;
    minimizerParams.A = 20;
    minimizerParams.a = 0.3;
    minimizerParams.SGDsamples = 2000;
    minimizerParams.n_walkers = 10;
    minimizerParams.thermalization = 100000;
    minimizerParams.n_cm = 1000;
    minimizerParams.n_c_SGD = 400;
    minimizerParams.alpha = arma::zeros(1, 1) + 0.5;
    minimizerParams.beta = arma::zeros(1, 1) + 0.5;


    //Seting values if not flagged default (controlled by Python)
    std::string def = "def";

    if (argc == n_args) {

        if (def.compare(argv[1]) != 0) outputParams.blocking_out = (bool)atoi(argv[1]);
        if (def.compare(argv[2]) != 0) outputParams.dist_out = (bool)atoi(argv[2]);
        if (def.compare(argv[3]) != 0) outputParams.outputSuffix = argv[3];
        if (def.compare(argv[4]) != 0) outputParams.outputPath = argv[4];

        if (def.compare(argv[5]) != 0) generalParams.n_p = atoi(argv[5]);
        if (def.compare(argv[6]) != 0) generalParams.dim = atoi(argv[6]);
        if (def.compare(argv[7]) != 0) generalParams.w = atof(argv[7]);

        if (def.compare(argv[8]) != 0) generalParams.random_seed = atoi(argv[8]);
        if (def.compare(argv[9]) != 0) generalParams.h = atof(argv[9]);
        if (def.compare(argv[10]) != 0) generalParams.D = atof(argv[10]);

        //        if (def.compare(argv[11]) != 0) parParams.parallel = (bool)atoi(argv[11]);

        if (def.compare(argv[12]) != 0) generalParams.doMIN = (bool)atoi(argv[12]);
        if (def.compare(argv[13]) != 0) generalParams.doVMC = (bool)atoi(argv[13]);
        if (def.compare(argv[14]) != 0) generalParams.doDMC = (bool)atoi(argv[14]);

        if (def.compare(argv[15]) != 0) generalParams.use_coulomb = (bool)atoi(argv[15]);
        if (def.compare(argv[16]) != 0) generalParams.use_jastrow = (bool)atoi(argv[16]);

        if (def.compare(argv[17]) != 0) generalParams.sampling = argv[17];
        //        if (def.compare(argv[18]) != 0) generalParams.kinetics_type = argv[18];
        if (def.compare(argv[19]) != 0) generalParams.system = argv[19];


        if (def.compare(argv[20]) != 0) vmcParams.n_c = atoi(argv[20]);
        if (def.compare(argv[21]) != 0) vmcParams.dt = atof(argv[21]);


        if (def.compare(argv[22]) != 0) dmcParams.dt = atof(argv[22]);
        if (def.compare(argv[23]) != 0) dmcParams.E_T = atof(argv[23]);
        if (def.compare(argv[24]) != 0) dmcParams.n_b = atoi(argv[24]);
        if (def.compare(argv[25]) != 0) dmcParams.n_w = atoi(argv[25]);
        if (def.compare(argv[26]) != 0) dmcParams.n_c = atoi(argv[26]);
        if (def.compare(argv[27]) != 0) dmcParams.therm = atoi(argv[27]);
        if (def.compare(argv[28]) != 0) dmcParams.dist_in = (bool)atoi(argv[28]);


        if (def.compare(argv[29]) != 0) minimizerParams.max_step = atof(argv[29]);
        if (def.compare(argv[30]) != 0) minimizerParams.f_max = atof(argv[30]);
        if (def.compare(argv[31]) != 0) minimizerParams.f_min = atof(argv[31]);
        if (def.compare(argv[32]) != 0) minimizerParams.omega = atof(argv[32]);
        if (def.compare(argv[33]) != 0) minimizerParams.A = atof(argv[33]);
        if (def.compare(argv[34]) != 0) minimizerParams.a = atof(argv[34]);
        if (def.compare(argv[35]) != 0) minimizerParams.SGDsamples = atoi(argv[35]);
        if (def.compare(argv[36]) != 0) minimizerParams.n_walkers = atoi(argv[36]);
        if (def.compare(argv[37]) != 0) minimizerParams.thermalization = atoi(argv[37]);
        if (def.compare(argv[38]) != 0) minimizerParams.n_cm = atoi(argv[38]);
        if (def.compare(argv[39]) != 0) minimizerParams.n_c_SGD = atoi(argv[39]);
        if (def.compare(argv[40]) != 0) minimizerParams.alpha = arma::zeros(1, 1) + atof(argv[40]);
        if (def.compare(argv[41]) != 0) minimizerParams.beta = arma::zeros(1, 1) + atof(argv[41]);

        if (def.compare(argv[42]) != 0) variationalParams.alpha = atof(argv[42]);
        if (def.compare(argv[43]) != 0) variationalParams.beta = atof(argv[43]);

    }


    //Seting timestep
    if (generalParams.sampling == "IS") {
        vmcParams.dt = 0.01;
    } else {
        vmcParams.dt = 0.5;
    }

    //Seting dmc inpath if dist out VMC is set
    if (outputParams.dist_out) {
        dmcParams.dist_in_path = outputParams.outputPath;
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

            if (dmcParams.n_w < parParams.n_nodes) {
                std::cout << "Unsufficient walkers for node structure." << std::endl;
            }

        }

        if (generalParams.doVMC && generalParams.doDMC) {
            if (dmcParams.dist_in && outputParams.dist_out) {
                if (dmcParams.n_w > vmcParams.n_c / (200)) {
                    std::cout << "Unsufficient VMC cycles to load dist in DMC." << std::endl;
                    std::cout << "For n_w=" << dmcParams.n_w << "the minimum is n_c=" << vmcParams.n_c / (200) << std::endl;
                }
            }
        }

        if (generalParams.doDMC) {

            if (dmcParams.dist_in && (!outputParams.dist_out)) {
                arma::mat r_test;
                std::stringstream s;
                s << dmcParams.dist_in_path << "walker_positions/dist_out0_0.arma";
                bool dataFound = r_test.load(s.str());
                if (!dataFound) {
                    std::cout << "No walker output data found in " << dmcParams.dist_in_path << "walker_positions/" << std::endl;
                    exit(1);
                }
                s.str(std::string());

                for (int i = 0; i < parParams.n_nodes; i++) {
                    s << dmcParams.dist_in_path << "walker_positions/dist_out" << i << "_0.arma";
                    dataFound = r_test.load(s.str());
                    if (!dataFound) {
                        std::cout << "No walker output data found for node " << i << ". Data generated on less nodes?" << std::endl;
                        exit(1);
                    }
                    s.str(std::string());
                }
                r_test.clear();

            }
        }
    }

    generalParams.random_seed -= parParams.node;
    minimizerParams.n_c_SGD /= parParams.n_nodes;
    vmcParams.n_c /= parParams.n_nodes;
    dmcParams.n_w /= parParams.n_nodes;

    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;
}