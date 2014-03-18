#include <QMC2.h>

#include <armadillo>

using namespace QMC2;
using namespace arma;

inline
void calcVirialPlot(int np, double w0, double w1, int Nw, ParParams & pp, double a0, double b0, std::string op);


int main(int argc, char** argv)
{

    (void) argc;
    ParParams pp;

    initMPI(pp, 0, NULL);

    std::string op = argv[1];
    op += "/";

    uint NW = 20;

    calcVirialPlot( 6, 0.01, 1, NW, pp,   0.9243,   0.5571, op);
    calcVirialPlot(12, 0.01, 1, NW, pp,   0.8756,     0.66, op);
    calcVirialPlot(20, 0.01, 1, NW, pp,   0.8361,   0.7332, op);
    calcVirialPlot(30, 0.01, 1, NW, pp,   0.8085,   0.7944, op);
    calcVirialPlot(42, 0.01, 1, NW, pp, 0.782778,  0.84400, op);
    calcVirialPlot(56, 0.01, 1, NW, pp,     0.76, 0.886972, op);


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
    mP.SGDsamples = 4000;

    mP.alpha(0) = a0;
    mP.beta(0) = b0;

    vmcP.n_c = pp.n_nodes*1E7;

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


    VMC vmc(gP, vmcP, sO, pp, 1, true);
    vmc.set_error_estimator(new SimpleVar(pp));

    ASGD asgd(&vmc, mP, pp, gP.runpath);

    CalculateAndSampleRadius sr;
    SampleRadiusSquared sr2;
    SampleRelativeDistance srij;

    vmc.add_subsample(&sr);
    vmc.add_subsample(&sr2);
    vmc.add_subsample(&srij);


    vec wList = linspace(w0, w1, Nw);

    mat results;
    std::stringstream name;

    if (pp.is_master)
    {
        name << gP.runpath << "/w_t_HO_COL_N" << gP.n_p << ".dat";
        results.set_size(Nw, 17);
    }

    for (int i = Nw-1; i >= 0; --i)
    {

        double w = wList(i);

        aHO.set_w(w);
        HO.set_w(w);

        bool init = i == Nw-1;
        asgd.minimize(init);
        vmc.run_method(init);

        double E = vmc.get_energy();
        double err_E = vmc.get_error();

        double T = vmc.kinetic_sampler.extract_mean();
        double vho = HO.pot_sampler.extract_mean();
        double vcol = COL.pot_sampler.extract_mean();
        double r = sr.extract_mean();
        double r2 = sr2.extract_mean();
        double rij = srij.extract_mean();


        double err_T = vmc.kinetic_sampler.extract_mean_error();
        double err_vho = HO.pot_sampler.extract_mean_error();
        double err_vcol = COL.pot_sampler.extract_mean_error();
        double err_r = sr.extract_mean_error();
        double err_r2 = sr2.extract_mean_error();
        double err_rij = srij.extract_mean_error();

        double alpha = sO.SP_basis->get_parameter(0);
        double beta  = sO.jastrow->get_parameter(0);

        if (pp.is_master)
        {
            results.row(i) = rowvec({w, alpha, beta, E, T, vho, vcol, r, r2, rij, err_E, err_T, err_vho, err_vcol, err_r, err_r2, err_rij});
            results.save(name.str(), raw_ascii);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

}
