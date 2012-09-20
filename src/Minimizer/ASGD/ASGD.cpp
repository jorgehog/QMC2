/* 
 * File:   ASGD.cpp
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#include "../../QMCheaders.h"

using namespace arma;

ASGD::ASGD(VMC* vmc,
        const rowvec & alpha,
        const rowvec & beta,
        int SGDsamples,
        int n_walkers,
        int n_c,
        int thermalization,
        int n_c_SGD,
        double max_step,
        double f_min,
        double f_max,
        double w,
        double a,
        double A)
: Minimizer(vmc, alpha, beta) {

    this->n_c = n_c;
    this->thermalization = thermalization;
    this->n_c_SGD = n_c_SGD;
    this->SGDsamples = SGDsamples;
    this->n_walkers = n_walkers;

    this->a = a;
    this->A = A;
    this->f_min = f_min;
    this->f_max = f_max;
    this->w = w;

    this->max_step = max_step;

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

    t_prev = 0;

}

double ASGD::f(double x) {
    return f_min + (f_max - f_min) / (1 - (f_max / f_min) * exp(-x / w));
}

void ASGD::get_variational_gradients(Walker* walker, double e_local) {

    for (int alpha = 0; alpha < Nspatial_params; alpha++) {

        double dalpha = vmc->get_orbitals_ptr()->get_variational_derivative(walker, alpha);

        gradient_local[alpha] += e_local * dalpha;
        gradient[alpha] += dalpha;

    }

    for (int beta = 0; beta < Njastrow_params; beta++) {

        double dbeta = vmc->get_jastrow_ptr()->get_variational_derivative(walker, beta);

        gradient_local[Nspatial_params + beta] += e_local*dbeta;
        gradient[Nspatial_params + beta] += dbeta;

    }

}

VMC* ASGD::minimize() {


    /* DEBUG */
    ofstream file;
    file.open("alpha.dat");
    DEBAG.open("sammenliknnBF_CF.dat");
    int debug1;
    double aGrad, bGrad, sumE;
    aGrad = bGrad = sumE = 0;
    debug1 = 1;
    //


    vmc->initialize();

    for (int cycle = 1; cycle <= thermalization + n_walkers * n_c; cycle++) {
        vmc->diffuse_walker(vmc->original_walker, vmc->trial_walker);

        if ((cycle > thermalization) && (cycle % n_c) == 0) {

            vmc->copy_walker(vmc->original_walker, walkers[k]);
            vmc->copy_walker(vmc->original_walker, trial_walkers[k]);
            k++;
        }

    }
    //D:100% match this far

    for (int sample = 1; sample <= SGDsamples; sample++) {

        E = 0;
        gradient = zeros(1, Nparams);
        gradient_local = zeros(1, Nparams);

        for (int k = 0; k < n_walkers; k++) {
            
            for (int cycle = 0; cycle < n_c_SGD; cycle++) {

                vmc->diffuse_walker(walkers[k], trial_walkers[k]);

                vmc->calculate_energy_necessities(walkers[k]);

                double e_local = vmc->calculate_local_energy(walkers[k]);
                E += e_local;

                get_variational_gradients(walkers[k], e_local);

            }
        }

        int scale = n_walkers*n_c_SGD;

        //        cout << E / scale << endl;
        E /= scale;

        gradient_tot = 2 * (gradient_local - gradient * E) / scale;

        double x = -dot(gradient_tot, gradient_old);
        double t_test = t_prev + f(x);

        t = t_test * (t_test > 0);

        //        output progress
        if ((sample % 100) == 0) {
            output("cycle:", sample);
            cout << gradient_tot << endl;
            cout << E << endl;
        }


        //        
        aGrad += gradient_tot[0];
        bGrad += gradient_tot[1];
        sumE += E;
        file << aGrad / debug1 << "\t" << bGrad / debug1 << "\t";
        //


        for (int param = 0; param < Nspatial_params; param++) {

            step = a / (t + A) * gradient_tot[param];
            if (step * step > max_step * max_step) {
                step *= max_step / fabs(step);
            }

            double alpha = vmc->get_orbitals_ptr()->get_parameter(param);
            file << abs(alpha - step) << "\t";
            vmc->get_orbitals_ptr()->set_parameter(abs(alpha - step), param);
        }

        for (int param = 0; param < Njastrow_params; param++) {

            step = a / (t + A) * gradient_tot[Nspatial_params + param];
            if (step * step > max_step * max_step) {
                step *= max_step / fabs(step);
            }

            double beta = vmc->get_jastrow_ptr()->get_parameter(param);
            file << abs(beta - step) << "\t";
            vmc->get_jastrow_ptr()->set_parameter(abs(beta - step), param);


        }

        file << sumE / debug1 << endl;
        debug1++;

        t_prev = t;
        gradient_old = gradient_tot;

    }

    output("Finished minimizing. Final parameters:", -1);
    file.close();
    DEBAG.close();
    vmc->accepted = 0;
    return vmc;
}