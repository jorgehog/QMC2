#include "ASGD.h"

#include "../../structs.h"

#include "../../Walker/Walker.h"
#include "../../Sampling/Sampling.h"
#include "../../System/System.h"
#include "../../Orbitals/Orbitals.h"
#include "../../Jastrow/Jastrow.h"
#include "../../OutputHandler/stdoutASGD/stdoutASGD.h"

#include <iomanip>


using namespace QMC2;

ASGD::ASGD(VMC* vmc, MinimizerParams & mP, const ParParams & pp, std::string path)
    : Minimizer(vmc, pp, mP.alpha, mP.beta) {
    using namespace arma;

    if (is_master) ASGDout = new stdoutASGD(this, path);
    
    this->n_c = mP.n_c;
    this->thermalization = mP.therm;
    this->n_c_SGD = mP.n_c_SGD;
    this->SGDsamples = mP.SGDsamples;
    this->n_walkers = mP.n_w;

    this->a = mP.a;
    this->A = mP.A;
    this->f_min = mP.f_min;
    this->f_max = mP.f_max;
    this->w = mP.omega;

    this->max_step = mP.max_step;

    walkers = new Walker*[n_walkers];
    trial_walkers = new Walker*[n_walkers];

    for (int i = 0; i < n_walkers; i++) {
        walkers[i] = new Walker(vmc->n_p, vmc->dim);
        trial_walkers[i] = new Walker(vmc->n_p, vmc->dim);
    }


    gradient = zeros(1, Nparams);
    gradient_local = zeros(1, Nparams);

    gradient_old = zeros(1, Nparams);
    gradient_tot = zeros(1, Nparams);

    t_prev.set_size(Nparams);
    t.set_size(Nparams);

    for (int i = 0; i < Nparams; ++i) {
        t_prev(i) = A;
    }
}

void ASGD::get_variational_gradients(Walker* walker, double e_local) {


    for (int alpha = 0; alpha < Nspatial_params; alpha++){
        double dalpha = vmc->get_orbitals_ptr()->get_variational_derivative(walker, alpha);

        gradient_local(alpha) += e_local * dalpha;
        gradient(alpha) += dalpha;
    }

    for (int beta = 0; beta < Njastrow_params; beta++) {

        double dbeta = vmc->get_jastrow_ptr()->get_variational_derivative(walker, beta);

        gradient_local(Nspatial_params + beta) += e_local*dbeta;
        gradient(Nspatial_params + beta) += dbeta;

    }

}

void ASGD::get_total_grad() {
    int scale = n_walkers * n_c_SGD;

#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &E, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    E /= n_nodes;
#endif
    E /= scale;

    gradient_tot = 2 * (gradient_local - gradient * E) / scale;


#ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, gradient_tot.memptr(), Nparams, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    gradient_tot /= n_nodes;
#endif


}

void ASGD::minimize(bool initialize) {

    vmc->get_sampling_ptr()->set_dt(vmc->dtOrig);

    if (initialize){
        initializeParameters();
        thermalize_walkers();
    } else {

        for (int i = 0; i < Nparams; ++i) {
            t_prev(i) = A;
        }

        gradient_tot.zeros();
        if (is_master) ASGDout->reset();
    }


    for (sample = 1; sample <= SGDsamples; sample++) {

        E = 0;
        gradient = arma::zeros(1, Nparams);
        gradient_local = arma::zeros(1, Nparams);

        for (int k = 0; k < n_walkers; k++) {
            vmc->get_sampling_ptr()->set_trial_states(walkers[k]);
            vmc->get_system_ptr()->initialize(walkers[k]);
            vmc->get_jastrow_ptr()->get_dJ_matrix(walkers[k]);
            vmc->get_sampling_ptr()->get_necessities(walkers[k]);

            for (int cycle = 0; cycle < n_c_SGD; cycle++) {

                vmc->diffuse_walker(walkers[k], trial_walkers[k]);

                vmc->calculate_energy_necessities(walkers[k]);
                double e_local = vmc->calculate_local_energy(walkers[k]);
                E += e_local;

                get_variational_gradients(walkers[k], e_local);

            }
        }

        get_total_grad();

        update_parameters();

        if (is_master) ASGDout->dump();

        output_cycle();

    }

    if (is_master) ASGDout->finalize();

    vmc->accepted = 0;
}

void ASGD::update_parameters() {

    for (int i = 0; i < Nparams; ++i) {

        double xi = -gradient_tot(i)*gradient_old(i);

        t(i) = t_prev(i) + f(xi);
        if (t(i) < 0) {
            t(i) = 0;
        }

    }


    for (int param = 0; param < Nspatial_params; param++) {

        step = a / (t(param) + A) * gradient_tot(param);
        if (fabs(step) > max_step) {
            step *= max_step / fabs(step);
        }

        double alpha = vmc->get_orbitals_ptr()->get_parameter(param);
        vmc->get_orbitals_ptr()->set_parameter(fabs(alpha - step), param);

    }

    for (int param = 0; param < Njastrow_params; param++) {

        step = a / (t(Nspatial_params + param) + A) * gradient_tot(Nspatial_params + param);
        if (step * step > max_step * max_step) {
            step *= max_step / fabs(step);
        }

        double beta = vmc->get_jastrow_ptr()->get_parameter(param);
        vmc->get_jastrow_ptr()->set_parameter(fabs(beta - step), param);

    }

    t_prev = t;
    gradient_old = gradient_tot;
}

void ASGD::output_cycle() {

    using namespace std;

    if (is_master) {

        if ((sample % 10) == 0) {



            if (Nspatial_params != 0) s << "a:";
            for (int alpha = 0; alpha < Nspatial_params; alpha++) {

                s << " ";
                s << setprecision(4) << fixed;
                s << vmc->get_orbitals_ptr()->get_parameter(alpha);

                s << " (";
                s << setprecision(4) << fixed << setfill(' ') << setw(7);
                s << gradient_tot(alpha) << ")";

            }

            if (Nspatial_params != 0) s << " | ";

            if (Njastrow_params != 0) s << "b:";
            for (int beta = 0; beta < Njastrow_params; beta++) {
                s << " ";
                s << setprecision(4) << fixed;
                s << vmc->get_jastrow_ptr()->get_parameter(beta);

                s << " (";
                s << setprecision(4) << fixed << setfill(' ') << setw(7);
                s << gradient_tot(Nspatial_params + beta) << ")";

            }

            s << " | E = ";
            s << setprecision(4) << fixed;
            s << E << " | ";


            s << setprecision(1) << fixed << setfill(' ') << setw(5);
            s << (double) sample / SGDsamples * 100 << "%";

            if (sample == SGDsamples) {
                cout << "\r" << s.str() << endl;
            } else {
                cout << "\r" << s.str();
                cout.flush();

#ifdef LINED_OUTPUT
                cout << endl;
#endif

            }

            s.str(string());
            s.clear();

        }
    }
}

void ASGD::thermalize_walkers() {
    int k = 0;
    vmc->get_sampling_ptr()->set_trial_pos(vmc->original_walker);
    vmc->copy_walker(vmc->original_walker, vmc->trial_walker);
    for (int cycle = 1; cycle <= thermalization + n_walkers * n_c; cycle++) {
        vmc->diffuse_walker(vmc->original_walker, vmc->trial_walker);

        if ((cycle > thermalization) && (cycle % n_c) == 0) {
            vmc->copy_walker(vmc->original_walker, walkers[k]);
            vmc->copy_walker(vmc->trial_walker, trial_walkers[k]);
            k++;
        }
    }
}
