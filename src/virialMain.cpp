#include "Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h"

#include "System/Fermions/Fermions.h"

#include "Potential/Harmonic_osc/Harmonic_osc.h"
#include "Potential/Coulomb/Coulomb.h"

#include "Jastrow/Pade_Jastrow/Pade_Jastrow.h"

#include "Sampling/Importance/Importance.h"

#include "structs.h"

#include "Minimizer/ASGD/ASGD.h"
#include "QMC/VMC/VMC.h"

#include "ErrorEstimator/SimpleVar/SimpleVar.h"


#include <armadillo>

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

    calcVirialPlot(2, 0.01, 1, 20, pp, 0.98831, 0.398664, op);
    calcVirialPlot(6, 0.01, 1, 20, pp, 0.9243, 0.5571, op);
    calcVirialPlot(12, 0.01, 1, 20, pp, 0.8756, 0.66, op);
    calcVirialPlot(20, 0.01, 1, 20, pp, 0.8361, 0.7332, op);
    calcVirialPlot(30, 0.01, 1, 20, pp, 0.8085, 0.7944, op);
    calcVirialPlot(42, 0.01, 1, 20, pp, 0.782778, 0.84400, op);
    calcVirialPlot(56, 0.01, 1, 20, pp, 0.76, 0.886972, op);


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
    mP.n_c_SGD = 25*pp.n_nodes;
    mP.SGDsamples = 1000;

    mP.alpha(0) = a0;
    mP.beta(0) = b0;
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

    vec wList = linspace(w0, w1, Nw);

    mat results;
    std::stringstream name;

    if (pp.is_master)
    {
        name << gP.runpath << "/w_t_HO_COL_N" << gP.n_p << ".dat";
        results.set_size(Nw, 4);
    }

    for (int i = Nw-1; i >= 0; --i)
    {

        double w = wList(i);

        aHO.set_w(w);
        HO.set_w(w);

        bool init = i == Nw-1;
        asgd.minimize(init);
        vmc.run_method(init);

        double T = vmc.kinetic_sampler.extract_mean();
        double vho = HO.pot_sampler.extract_mean();
        double vcol = COL.pot_sampler.extract_mean();

        if (pp.is_master)
        {
            results.row(i) = rowvec({w, T, vho, vcol});
            results.save(name.str(), raw_ascii);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

}