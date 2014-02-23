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
void calcVirialPlot(int np, double w0, double w1, int Nw, ParParams & pp, double a0, double b0);


int main()
{

    ParParams pp;

    initMPI(pp, 0, NULL);

    calcVirialPlot(2, 0.01, 1, 20, pp, 0.98831, 0.398664);
    calcVirialPlot(6, 0.01, 1, 20, pp, 0.9243, 0.5571);
    calcVirialPlot(12, 0.01, 1, 20, pp, 0.8756, 0.66);
    calcVirialPlot(20, 0.01, 1, 20, pp, 0.8361, 0.7332);
    calcVirialPlot(30, 0.01, 1, 20, pp, 0.8085, 0.7944);
    calcVirialPlot(42, 0.01, 1, 20, pp, 0.782778, 0.84400);
    calcVirialPlot(56, 0.01, 1, 20, pp, 0.76, 0.886972);


    MPI_Finalize();

    return 0;

}

inline
void calcVirialPlot(int np, double w0, double w1, int Nw, ParParams & pp, double a0, double b0)
{

    GeneralParams gP;
    SystemObjects sO;
    VariationalParams vP;
    MinimizerParams mP;
    VMCparams vmcP;
    DMCparams dmcP;

    gP.n_p = np;
    mP.n_c_SGD = 100;
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

    if (pp.is_master)
    {
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
        }

    }

    if (pp.is_master)
    {

        std::stringstream name;
        name << gP.runpath << "/w_t_HO_COL_N" << gP.n_p << ".dat";

        results.save(name.str(), raw_ascii);
    }

    MPI_Barrier(MPI_COMM_WORLD);

}
