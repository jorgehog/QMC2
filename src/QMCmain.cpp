/* 
 * File:   QMCmain.cpp
 * Author: jorgehog
 *
 * Created on 13. april 2012, 17:04
 */

//#include "mpi.h"
#include "QMCheaders.h"

#include <sys/time.h>

/*
 * 
 */

void parseCML(char**,
        VMCparams &,
        DMCparams &,
        VariationalParams &,
        GeneralParams &,
        MinimizerParams &,
        OutputParams &,
        SystemObjects &
        );

int main(int argc, char** argv) {
    using namespace std;

    struct VMCparams vmcParams;
    struct DMCparams dmcParams;
    struct VariationalParams variationalParams;
    struct GeneralParams generalParams;
    struct MinimizerParams minimizerParams;
    struct OutputParams outputParams;
    struct SystemObjects systemObjects;
    
    int n_p, dim, n_c, numprocs, my_rank;
    double alpha, beta, w, dt, h, cumul_e, cumul_e2, e, e2, E_T;
    long random_seed;

    //    parseCML(argv, vmcParams,
    //            dmcParams,
    //            variationalParams,
    //            generalParams,
    //            minimizerParams,
    //            outputParams,
    //            systemObjects);

    //initializing MPI
    //    MPI_Init(&argc, &argv);
    //    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    random_seed = -1.234;
    //    random_seed = -time(NULL);
    //    random_seed = -time(NULL) - my_rank;

    n_p = 6;
    dim = 2;
    w = 1;

    n_c = 1E6;
    bool doMin = false;
    bool doVmc = false;
    bool doDmc = true;

    string system = "QDots";
    string sampling = "BF";
    string kinetics_type = "CF";

    bool dist_out = true;
    bool blocking_out = false;

    bool dist_in = true;

    bool use_jastrow = true;
    bool use_coulomb = true;


    initVMC(n_p, dim, w, dt, system, sampling, alpha, beta);
    //cout << alpha << " " << beta << endl;
    if ((use_jastrow == false) && (use_coulomb == false)) {
        alpha = 1;
    }


    Kinetics* kinetics;
    Orbitals* SP_basis;
    Potential* onebody_pot;
    System* SYSTEM;
    Sampling* sample_method;
    Jastrow* jastrow;

    if (kinetics_type == "Num") {
        h = 0.0001;
        kinetics = new Numerical(n_p, dim, h);
    } else if (kinetics_type == "CF") {
        kinetics = new Closed_form(n_p, dim);
    } else {
        cout << "unknown kinetics" << endl;
        exit(1);
    }

    if (sampling == "IS") {
        sample_method = new Importance(n_p, dim, dt, random_seed);
    } else if (sampling == "BF") {
        sample_method = new Brute_Force(n_p, dim, dt, random_seed);
    } else {
        cout << "unknown sampling method" << endl;
        exit(1);
    }

    if (use_jastrow) {
        jastrow = new Pade_Jastrow(n_p, dim, beta);
    } else {
        jastrow = new No_Jastrow();
    }

    if (system == "QDots") {
        SP_basis = new oscillator_basis(n_p, dim, alpha, w);

        onebody_pot = new Harmonic_osc(n_p, dim, w);


        SYSTEM = new Fermions(n_p, dim, SP_basis);
        SYSTEM->add_potential(onebody_pot);


    } else {
        cout << "unknown system" << endl;
        exit(1);
    }

    if (use_coulomb) {
        SYSTEM->add_potential(new Coulomb(n_p, dim));
    }



    if (doVmc || doMin) {

        VMC* vmc = new VMC(n_p, dim, n_c,
                sample_method,
                SYSTEM,
                kinetics,
                jastrow);
        if (dist_out) {

            OutputHandler* dist = new Distribution("dist_out");
            vmc->add_output(dist);

        }

        if (blocking_out) {

            OutputHandler* blocking = new BlockingData("blocking_out");
            vmc->add_output(blocking);

        }

        if (doMin) {
            double max_step = 0.1;
            double f_max = 1.0;
            double f_min = -0.5;
            double omega = 0.8;
            double A = 20;
            double a = 0.3;
            int SGDsamples = 10000;
            int n_walkers = 10;
            int thermalization = 100000;
            int n_cm = 1000;
            int n_c_SGD = 100;
            rowvec alpha = zeros(1, 1) + 0.5;
            rowvec beta = zeros(1, 1) + 0.5;



            //            Minimizer * minimizer = new ASGD(vmc,
            //                    alpha,
            //                    beta,
            //                    SGDsamples,
            //                    n_walkers,
            //                    n_cm,
            //                    thermalization,
            //                    n_c_SGD,
            //                    max_step,
            //                    f_min,
            //                    f_max,
            //                    omega,
            //                    a,
            //                    A);
            //
            //            vmc = minimizer->minimize();
        }

        if (doVmc) {
            vmc->run_method();
        }


        //    cumul_e = cumul_e2 = 0;
        //    e = vmc->get_energy();
        //    e2 = vmc->get_e2();
        //    
        //    MPI_Reduce(&e, &cumul_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //    MPI_Reduce(&e2, &cumul_e2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //
        //    if (my_rank == 0) {
        //        vmc->set_e(e);
        //        vmc->set_e2(e2);
        //        
        //        vmc->output();
        //    }

        if (doDmc) {
            E_T = vmc->get_energy();
        }
    }
    if (doDmc) {

        n_c = 1000;
        int therm = 1000;
        dt = 0.001;
        int n_w = 1000;
        int n_b = 200;

        sample_method->set_dt(dt);

        DMC* dmc = new DMC(n_p, dim, n_w, n_c, n_b, therm, E_T,
                sample_method,
                SYSTEM,
                kinetics,
                jastrow,
                dist_in);

        OutputHandler* DMCout = new stdoutDMC();

        dmc->add_output(DMCout);

        dmc->run_method();
        cout << "DMC FIN." << endl;
    }



    //    MPI_Finalize();
    //    ofstream paramSpaceOut;
    //    paramSpaceOut.open("alphaVsE.dat");
    //    VMC* vmc;
    //    for (double alpha = 0.5; alpha <= 1.5; alpha+=0.1){
    //        SP_basis = new oscillator_basis(n_p, dim, alpha, w);
    //        SYSTEM = new Fermions(n_p, dim, SP_basis);
    //        SYSTEM->add_potential(onebody_pot);
    //        vmc = new VMC(n_p, dim, n_c,
    //                sample_method,
    //                SYSTEM,
    //                kinetics,
    //                jastrow);
    //        cout << alpha << endl;
    //        vmc->run_method();
    //        paramSpaceOut << alpha << "\t" << vmc->get_energy() << endl;
    //    }
    //    paramSpaceOut.close();
    return 0;
}

void parseCML(char** argv, VMCparams & vmcParams,
        DMCparams & dmcParams,
        VariationalParams & variationalParams,
        GeneralParams & generalParams,
        MinimizerParams & minimizerParams,
        OutputParams & outputParams,
        SystemObjects & systemObjects) {

    //Default values:
    generalParams.n_p = 2;
    generalParams.dim = 2;
    generalParams.w = 1;

    generalParams.random_seed = -time(NULL);
    generalParams.h = 0.0001;

    generalParams.parallell = false;
    generalParams.numprocs = 1;
    generalParams.myrank = 0;

    generalParams.doMIN = true;
    generalParams.doVMC = true;
    generalParams.doDMC = true;

    generalParams.use_coulomb = true;
    generalParams.use_jastrow = true;

    generalParams.sampling = "IS";
    generalParams.kinetics_type = "CF";
    generalParams.system = "QDots";


    //    vmcParams.dt = 0.5;
    vmcParams.n_c = 10E6;
    vmcParams.therm = 10E5;


    dmcParams.dt = 0.001;
    dmcParams.E_T = 0;
    dmcParams.n_b = 100;
    dmcParams.n_w = 1000;
    dmcParams.n_c = 1000;
    dmcParams.therm = 1000;
    dmcParams.dist_in = false;

    minimizerParams.max_step = 0.1;
    minimizerParams.f_max = 1.0;
    minimizerParams.f_min = -0.5;
    minimizerParams.omega = 0.8;
    minimizerParams.A = 20;
    minimizerParams.a = 0.3;
    minimizerParams.SGDsamples = 10000;
    minimizerParams.n_walkers = 10;
    minimizerParams.thermalization = 100000;
    minimizerParams.n_cm = 1000;
    minimizerParams.n_c_SGD = 100;
    minimizerParams.alpha = zeros(1, 1) + 0.5;
    minimizerParams.beta = zeros(1, 1) + 0.5;


    outputParams.blocking_out = false;
    outputParams.dist_out = false;
    outputParams.outputSuffix = "";
    outputParams.outputPath = "./";





    if (generalParams.kinetics_type == "Num") {
        systemObjects.kinetics = new Numerical(generalParams);
    } else if (generalParams.kinetics_type == "CF") {
        systemObjects.kinetics = new Closed_form(generalParams);
    } else {
        cout << "unknown kinetics" << endl;
        exit(1);
    }

    if (generalParams.sampling == "IS") {
        systemObjects.sample_method = new Importance(generalParams);
    } else if (generalParams.sampling == "BF") {
        systemObjects.sample_method = new Brute_Force(generalParams);
    } else {
        cout << "unknown sampling method" << endl;
        exit(1);
    }


    if (generalParams.use_jastrow) {
        systemObjects.jastrow = new Pade_Jastrow(generalParams, variationalParams);
    } else {
        systemObjects.jastrow = new No_Jastrow();
    }

    if (generalParams.system == "QDots") {
        systemObjects.SP_basis = new oscillator_basis(generalParams, variationalParams);

        systemObjects.onebody_pot = new Harmonic_osc(generalParams);

        systemObjects.SYSTEM = new Fermions(generalParams, systemObjects.SP_basis);
        systemObjects.SYSTEM->add_potential(systemObjects.onebody_pot);


    } else {
        cout << "unknown system" << endl;
        exit(1);
    }

    if (generalParams.use_coulomb) {
        systemObjects.SYSTEM->add_potential(new Coulomb(generalParams));
    }



    if (generalParams.doVMC || generalParams.doMIN) {

        VMC* vmc = new VMC(generalParams, vmcParams, systemObjects);
        systemObjects.sample_method->set_dt(vmcParams.dt);

        if (outputParams.dist_out) {
            string distname = "dist_out" + outputParams.outputSuffix;
            OutputHandler* dist = new Distribution(distname, outputParams.outputPath,
                    generalParams.parallell, generalParams.myrank, generalParams.numprocs);
            vmc->add_output(dist);

        }

        if (outputParams.blocking_out) {
            string blockname = "blocking_out" + outputParams.outputSuffix;
            OutputHandler* blocking = new BlockingData(blockname, outputParams.outputPath,
                    generalParams.parallell, generalParams.myrank, generalParams.numprocs);
            vmc->add_output(blocking);

        }

        if (generalParams.doMIN) {

            Minimizer * minimizer = new ASGD(vmc, minimizerParams);

            vmc = minimizer->minimize();
        }

        if (generalParams.doVMC) {
            vmc->run_method();
        }


        //    cumul_e = cumul_e2 = 0;
        //    e = vmc->get_energy();
        //    e2 = vmc->get_e2();
        //    
        //    MPI_Reduce(&e, &cumul_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //    MPI_Reduce(&e2, &cumul_e2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        //
        //    if (my_rank == 0) {
        //        vmc->set_e(e);
        //        vmc->set_e2(e2);
        //        
        //        vmc->output();
        //    }

    }

    if (generalParams.doDMC) {


        systemObjects.sample_method->set_dt(dmcParams.dt);

        DMC* dmc = new DMC(generalParams, dmcParams, systemObjects);

        string dmcOutname = (string) "DMC_out" + outputParams.outputSuffix;
        OutputHandler* DMCout = new stdoutDMC(dmcOutname, outputParams.outputPath,
                generalParams.parallell, generalParams.myrank, generalParams.numprocs);

        dmc->add_output(DMCout);

        dmc->run_method();

        cout << "DMC FIN." << endl;
    }
}