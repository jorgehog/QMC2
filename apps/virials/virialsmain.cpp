#include <QMC2.h>

#include <map>

#include <armadillo>

using namespace QMC2;
using namespace arma;
using namespace std;

inline
void calcVirialPlot(int np, double w0, double w1, int Nw, ParParams & pp, double a0, double b0, std::string op);


int main(int argc, char** argv)
{

    ParParams pp;

    initMPI(pp, argc, argv);

    if (argc != 3)
    {
        cerr << "Supply nparticles as first and output path as second CML arg." << endl;
        exit(1);
    }

    int n_p = atoi(argv[1]);
    string op = argv[2];
    op += "/";

    uint NW = 5;

    Parameters p;

    Parameters::system_maptype systemParameters = p.fetch_system("QDots");

    const double alpha = systemParameters[n_p][1.0][0];
    const double beta = systemParameters[n_p][1.0][1];

    calcVirialPlot(n_p, 0.01, 1, NW, pp, alpha, beta, op);
//    calcVirialPlot(12, 0.01, 1, NW, pp, 0.8756, 0.66, op);
//    calcVirialPlot(20, 0.01, 1, NW, pp, 0.8361, 0.7332, op);
//    calcVirialPlot(30, 0.01, 1, NW, pp, 0.8085, 0.7944, op);
//    calcVirialPlot(42, 0.01, 1, NW, pp, 0.782778, 0.84400, op);
//    calcVirialPlot(56, 0.01, 1, NW, pp, 0.76, 0.886972, op);


    MPI_Finalize();

    return 0;

}

inline
void calcVirialPlot(int np, double w0, double w1, int Nw, ParParams & pp, double a0, double b0, std::string op)
{

    GeneralParams gP;
    SystemObjects sO;
    VariationalParams vP;
    MinimizerParams mP;
    VMCparams vmcP;
    DMCparams dmcP;

    gP.runpath = op;

    gP.n_p = np;
    mP.n_c_SGD = 100*pp.n_nodes;
    mP.SGDsamples = 2000;

    mP.alpha(0) = a0;
    mP.beta(0) = b0;

    vmcP.n_c = pp.n_nodes*1E6;

    dmcP.n_w = pp.n_nodes*250;

    scaleWithProcs(pp, gP, mP, vmcP, dmcP);

    AlphaHarmonicOscillator aHO(gP, vP);
    sO.SP_basis = &aHO;

    sO.system = new Fermions(gP, sO.SP_basis);

    Harmonic_osc HO(gP);
    Coulomb COL(gP);

    sO.system->add_potential(&HO);
    sO.system->add_potential(&COL);

    sO.jastrow = new Pade_Jastrow(gP, vP);

    sO.sample_method = new Importance(gP);


    VMC vmc(gP, vmcP, sO, pp, dmcP.n_w, true);
    vmc.set_error_estimator(new Blocking(vmcP.n_c, pp, "blocking_VMC_out", gP.runpath));

    ASGD asgd(&vmc, mP, pp, gP.runpath);

    DMC dmc(gP, dmcP, sO, pp, true);
    dmc.set_error_estimator(new Blocking(dmcP.n_c, pp, "blocking_DMC_out", gP.runpath));

    CalculateAndSampleRadius sr_vmc;
    SampleRadiusSquared sr2_vmc;
    SampleRelativeDistance srij_vmc;

    vmc.add_subsample(&sr_vmc);
    vmc.add_subsample(&sr2_vmc);
    vmc.add_subsample(&srij_vmc);

    CalculateAndSampleRadius sr_dmc;
    SampleRadiusSquared sr2_dmc;
    SampleRelativeDistance srij_dmc;

    dmc.add_subsample(&sr_dmc);
    dmc.add_subsample(&sr2_dmc);
    dmc.add_subsample(&srij_dmc);



    vec wList = linspace(w0, w1, Nw);

    mat results;
    std::stringstream name;

    if (pp.is_master)
    {
        name << gP.runpath << "/w_t_HO_COL_N" << gP.n_p << ".dat";
        results.set_size(Nw, 31);
    }

    results.zeros();


    for (int i = Nw-1; i >= 0; --i)
    {

        double w = wList(i);

        stringstream s;
        s << "N" << gP.n_p << "w" << w;

        OutputHandler::suffix = s.str();

        aHO.set_w(w);
        HO.set_w(w);

        bool init = i == Nw-1;


        //MINIMIZING
        asgd.minimize(init);

        double alpha = sO.SP_basis->get_parameter(0);
        double beta  = sO.jastrow->get_parameter(0);

        //VMC
        vmc.run_method(init);

        double E_vmc = vmc.get_energy();
        double err_E_vmc = vmc.get_error();

        double T_vmc = vmc.kinetic_sampler.extract_mean();
        double vho_vmc = HO.pot_sampler.extract_mean();
        double vcol_vmc = COL.pot_sampler.extract_mean();
        double r_vmc = sr_vmc.extract_mean();
        double r2_vmc = sr2_vmc.extract_mean();
        double rij_vmc = srij_vmc.extract_mean();

        double err_T_vmc = vmc.kinetic_sampler.extract_mean_error();
        double err_vho_vmc = HO.pot_sampler.extract_mean_error();
        double err_vcol_vmc = COL.pot_sampler.extract_mean_error();
        double err_r_vmc = sr_vmc.extract_mean_error();
        double err_r2_vmc = sr2_vmc.extract_mean_error();
        double err_rij_vmc = srij_vmc.extract_mean_error();

        //DMC
        if (init)
        {
            dmc.initFromVMC(&vmc);
            vmc.do_store_walkers = false;
        }

        dmc.run_method(init);

        double E_dmc = dmc.get_energy();
        double err_E_dmc = dmc.get_error();

        double T_dmc = dmc.kinetic_sampler.extract_mean_of_means();
        double vho_dmc = HO.pot_sampler.extract_mean_of_means();
        double vcol_dmc = COL.pot_sampler.extract_mean_of_means();
        double r_dmc = sr_dmc.extract_mean_of_means();
        double r2_dmc = sr2_dmc.extract_mean_of_means();
        double rij_dmc = srij_dmc.extract_mean_of_means();

        double err_T_dmc = dmc.kinetic_sampler.extract_mean_of_means_error();
        double err_vho_dmc = HO.pot_sampler.extract_mean_of_means_error();
        double err_vcol_dmc = COL.pot_sampler.extract_mean_of_means_error();
        double err_r_dmc = sr_dmc.extract_mean_of_means_error();
        double err_r2_dmc = sr2_dmc.extract_mean_of_means_error();
        double err_rij_dmc = srij_dmc.extract_mean_of_means_error();


        //dump data

        if (pp.is_master)
        {
            results.row(i) = rowvec({w,
                                     alpha,
                                     beta,
                                     E_vmc,
                                     T_vmc,
                                     vho_vmc,
                                     vcol_vmc,
                                     r_vmc,
                                     r2_vmc,
                                     rij_vmc,
                                     err_E_vmc,
                                     err_T_vmc,
                                     err_vho_vmc,
                                     err_vcol_vmc,
                                     err_r_vmc,
                                     err_r2_vmc,
                                     err_rij_vmc,
                                     E_dmc,
                                     T_dmc,
                                     vho_dmc,
                                     vcol_dmc,
                                     r_dmc,
                                     r2_dmc,
                                     rij_dmc,
                                     err_E_dmc,
                                     err_T_dmc,
                                     err_vho_dmc,
                                     err_vcol_dmc,
                                     err_r_dmc,
                                     err_r2_dmc,
                                     err_rij_dmc});

            results.save(name.str(), raw_ascii);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

}
