#include <QMC2.h>

#include <map>

#include <armadillo>

#include <iomanip>

using namespace QMC2;
using namespace arma;
using namespace std;

inline
void calcVirialPlot(int np, double w0, double w1, int Nw, ParParams & pp, double a0, double b0, std::string op);

class EnergySampler : public Sampler
{
public:
    EnergySampler() : Sampler("Energy") {}

    // Sampler interface
public:
    void push_values(const Walker *walker)
    {
        push_value(walker->E);
    }
};


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

    uint NW = 20;

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
    mP.SGDsamples = 500;

    mP.alpha(0) = a0;
    mP.beta(0) = b0;

    vmcP.n_c = pp.n_nodes*1E6;

    dmcP.n_w = pp.n_nodes*250;
    dmcP.therm = 500;


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
    SampleRelativeDistance srij;

    vmc.add_subsample(&sr);
    vmc.add_subsample(&srij);

    dmc.add_subsample(&sr);
    dmc.add_subsample(&srij);

    vec wList = linspace(w0, w1, Nw);

    mat results;
    std::stringstream name;

    if (pp.is_master)
    {
        name << gP.runpath << "/w_t_HO_COL_N" << gP.n_p << ".dat";
        results.set_size(Nw, 15);
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

        if (pp.is_master)
        {
            cout << "\nRunning w = " << w  << " (" << Nw - i << "/" << Nw << ")" << endl;
        }

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
        double rij_vmc = srij.result();

        double err_T_vmc = vmc.kinetic_sampler->error();
        double err_vho_vmc = HO.pot_sampler->error();
        double err_vcol_vmc = COL.pot_sampler->error();
        double err_r_vmc = sr.error();
        double err_rij_vmc = srij.error();

        if (pp.is_master)
        {
            cout << "<E>    " << setprecision(8) << fixed << E_vmc    << " " << err_E_vmc << endl;
            cout << "<T>    " << setprecision(8) << fixed << T_vmc    << " " << err_T_vmc << endl;
            cout << "<vho>  " << setprecision(8) << fixed << vho_vmc  << " " << err_vho_vmc << endl;
            cout << "<vcol> " << setprecision(8) << fixed << vcol_vmc << " " << err_vcol_vmc << endl;
            cout << "<r>    " << setprecision(8) << fixed << r_vmc    << " " << err_r_vmc << endl;
            cout << "<rij>  " << setprecision(8) << fixed << rij_vmc  << " " << err_rij_vmc << endl;
        }

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
        double rij_dmc = srij.result();

        double err_T_dmc = dmc.kinetic_sampler->error();
        double err_vho_dmc = HO.pot_sampler->error();
        double err_vcol_dmc = COL.pot_sampler->error();
        double err_r_dmc = sr.error();
        double err_rij_dmc = srij.error();

        if (pp.is_master)
        {
            cout << "<E>    " << setprecision(8) << fixed << E_dmc    << " " << err_E_dmc   << endl;
            cout << "<T>    " << setprecision(8) << fixed << T_dmc    << " " << err_T_dmc   << endl;
            cout << "<vho>  " << setprecision(8) << fixed << vho_dmc  << " " << err_vho_dmc << endl;
            cout << "<vcol> " << setprecision(8) << fixed << vcol_dmc << " " << err_vcol_dmc<< endl;
            cout << "<r>    " << setprecision(8) << fixed << r_dmc    << " " << err_r_dmc   << endl;
            cout << "<rij>  " << setprecision(8) << fixed << rij_dmc  << " " << err_rij_dmc << endl;
        }

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
                                     rij_vmc,
                                     E_dmc,
                                     T_dmc,
                                     vho_dmc,
                                     vcol_dmc,
                                     r_dmc,
                                     rij_dmc});

            results.save(name.str(), raw_ascii);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

}
