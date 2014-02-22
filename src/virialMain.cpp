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

int main()
{

    GeneralParams gP;
    SystemObjects sO;
    VariationalParams vP;
    MinimizerParams mP;
    ParParams pp;
    VMCparams vmcP;
    DMCparams dmcP;

    initMPI(pp, 0, NULL);
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


    uint N = 20;
    vec wList = linspace(0.01, 1, N);

    mat results;

    if (pp.is_master)
    {
        results.set_size(N, 4);
    }

    for (int i = N-1; i >= 0; --i)
    {

        double w = wList(i);

        aHO.set_w(w);
        HO.set_w(w);

        bool init = (uint)i == N-1;
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

    MPI_Finalize();

    return 0;

}
