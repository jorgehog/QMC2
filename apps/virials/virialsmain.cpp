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

    ErrorEstimator::default_path = op;
    Sampler::standardErrorEstimator = Sampler::BLOCKING;

    stringstream ss;
    ss << "N" << gP.n_p << "w" << w1;
    OutputHandler::suffix = ss.str();

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
    vmc.set_error_estimator(new Blocking(vmcP.n_c, pp, "blocking_out", gP.runpath));

    ASGD asgd(&vmc, mP, pp, gP.runpath);

    DMC dmc(gP, dmcP, sO, pp, true);
    dmc.set_error_estimator(new Blocking(dmcP.n_c, pp, "blocking_out", gP.runpath));

    CalculateAndSampleRadius sr;
    SampleRadiusSquared sr2;
    SampleRelativeDistance srij;

    vmc.add_subsample(&sr);
    vmc.add_subsample(&sr2);
    vmc.add_subsample(&srij);

    dmc.add_subsample(&sr);
    dmc.add_subsample(&sr2);
    dmc.add_subsample(&srij);



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

        double T_vmc = vmc.kinetic_sampler->result();
        double vho_vmc = HO.pot_sampler->result();
        double vcol_vmc = COL.pot_sampler->result();
        double r_vmc = sr.result();
        double r2_vmc = sr2.result();
        double rij_vmc = srij.result();

        double err_T_vmc = vmc.kinetic_sampler->error();
        double err_vho_vmc = HO.pot_sampler->error();
        double err_vcol_vmc = COL.pot_sampler->error();
        double err_r_vmc = sr.error();
        double err_r2_vmc = sr2.error();
        double err_rij_vmc = srij.error();

        //DMC
        if (init)
        {
            dmc.initFromVMC(&vmc);
            vmc.do_store_walkers = false;
        }

        dmc.run_method(init);

        double E_dmc = dmc.get_energy();
        double err_E_dmc = dmc.get_error();

        double T_dmc = dmc.kinetic_sampler->result();
        double vho_dmc = HO.pot_sampler->result();
        double vcol_dmc = COL.pot_sampler->result();
        double r_dmc = sr.result();
        double r2_dmc = sr2.result();
        double rij_dmc = srij.result();

        double err_T_dmc = dmc.kinetic_sampler->result();
        double err_vho_dmc = HO.pot_sampler->result();
        double err_vcol_dmc = COL.pot_sampler->result();
        double err_r_dmc = sr.result();
        double err_r2_dmc = sr2.result();
        double err_rij_dmc = srij.result();


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
