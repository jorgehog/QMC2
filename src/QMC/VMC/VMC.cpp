#include "VMC.h"

#include <iomanip>

#include "../../structs.h"

#include "../../Walker/Walker.h"
#include "../../Sampling/Sampling.h"
#include "../../ErrorEstimator/ErrorEstimator.h"
#include "../../OutputHandler/Distribution/Distribution.h"
#include "../../Orbitals/Orbitals.h"
#include "../../System/System.h"


using namespace QMC2;


VMC::VMC(GeneralParams & gP, VMCparams & vP, SystemObjects & sO, ParParams & pp, int n_w, bool silent)
: QMC(gP, vP.n_c, sO, pp, silent, vP.dt, n_w) {

    std::stringstream name;
    name << sO.SP_basis->getName() << "_vmc";
    distribution = new Distribution(pp, gP.runpath, name.str(), silent);
    kinetic_sampler = new Sampler("kinetic");

    pop_tresh = n_c / n_w;
    offset = n_c - n_w*pop_tresh;
    dist_tresh = 25;
    last_walker = 0;

    dist_size = n_c * n_p / dist_tresh;
    dist = arma::zeros(dist_size, dim);

    output_tresh = n_c / 100;

    if (output_tresh == 0) {
        output_tresh = 1;
    }

    original_walker = new Walker(n_p, dim);

    vmc_E = 0;
    thermalization = n_c / 10 * (n_c < 1e6) + 1e5 * (n_c >= 1e6);

}

VMC::~VMC() {
    delete trial_walker;
    delete [] original_walkers;
}

void VMC::set_trial_positions() {
    sampling->set_trial_pos(original_walker);
}

void VMC::save_distribution() {
    using namespace arma;

    if (cycle % dist_tresh == 0) {
        dist(span(last_inserted, last_inserted + n_p - 1), span()) = original_walker->r;
        last_inserted += n_p;
    }
}

void VMC::store_walkers() {

    if (!do_store_walkers)
    {
        return;
    }

    if ((cycle > offset) && (cycle % pop_tresh == 0)) {
        copy_walker(original_walker, original_walkers[last_walker]);
        original_walkers[last_walker]->set_E(local_E);
        last_walker++;
    }

}


void VMC::run_method(bool initialize) {

    initializeRun("VMC");

    sampling->set_dt(dtOrig);
    m_currentlyRunningMethod = "VMC";

    kinetic_sampler->getMeanErrorEstimator()->setNumberOfCycles(n_c);
    system->setMeanErrorEstimatorNumberOfCycles(n_c);

    for (Sampler * sampler_method : samplers) {
        sampler_method->getMeanErrorEstimator()->setNumberOfCycles(n_c);
    }

    if (initialize) {

        set_trial_positions();

        copy_walker(original_walker, trial_walker);

        for (cycle = 1; cycle <= thermalization; cycle++) {
            diffuse_walker(original_walker, trial_walker);
        }

    } else {
        reset_all();
    }



    for (cycle = 1; cycle <= n_c; cycle++) {

        diffuse_walker(original_walker, trial_walker);

        calculate_energy_necessities(original_walker);

        local_E = calculate_local_energy(original_walker);

        vmc_E += local_E;

        update_samplers(original_walker, 1);

        output();
        update_subsamples();
        save_distribution();
        error_estimator->update_data(local_E);
        store_walkers();

    }

    node_comm();
    output();

    get_accepted_ratio();

    finalize_distribution();
    estimate_error();

    finalize();

    clean();
}

void VMC::output() {

    using namespace std;

    stringstream s;

    if (is_master) {

        if (cycle == n_c + 1) {

            s <<  "VMC energy: ";
            s << setprecision(6) << fixed << vmc_E;
            s << "  ";
            s << setprecision(1) << fixed << setfill(' ') << setw(5);
            s << 100.0 << "%";

            cout << "\r" << s.str() << endl;

            s.str(string());
            s.clear();

        } else if (cycle % output_tresh == 0) {

            s << "VMC energy: ";
            s << setprecision(6) << fixed << vmc_E / cycle;
            s << "  ";
            s << setprecision(1) << fixed << setfill(' ') << setw(5);
            s << (double) cycle / n_c * 100 << "%";

            cout << "\r" << s.str();
            cout.flush();

#ifdef LINED_OUTPUT
            cout << endl;
#endif
        }

    }

}

void VMC::node_comm() {
#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &vmc_E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    scale_values();

}

void VMC::reset_all()
{

    last_walker = 0;
    vmc_E = 0;

    QMC::reset_all();
}


