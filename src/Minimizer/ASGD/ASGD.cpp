/* 
 * File:   ASGD.cpp
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#include "../../QMCheaders.h"

using namespace arma;

ASGD::ASGD(VMC* vmc, int SGDsamples, int n_walkers, int n_c, int n_c_SGD, double f_min, double f_max, double w,
        double a, double A, int NSP, int NJP)
: Minimizer(vmc, NSP, NJP) {

    this->n_c = n_c;
    this->n_c_SGD = n_c_SGD;
    this->SGDsamples = SGDsamples;
    this->n_walkers = n_walkers;
    
    this->a = a;
    this->A = A;
    this->f_min = f_min;
    this->f_max = f_max;
    this->w = w;

    walkers = new Walker*[n_walkers];
    for (int i = 0; i < n_walkers; i++) {
        walkers[i] = new Walker(vmc->n_p, vmc->dim, false);
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

void ASGD::get_variational_gradients(double e_local) {

    for (int alpha = 0; alpha < Nspatial_params; alpha++) {
        double dalpha = vmc->get_orbitals_ptr()->get_variational_derivative(vmc->original_walker, alpha);

        gradient_local[alpha] += e_local * dalpha;
        gradient[alpha] += dalpha;

    }



    for (int beta = 0; beta < Njastrow_params; beta++) {

        double dbeta = vmc->get_jastrow_ptr()->get_variational_derivative(vmc->original_walker, beta);

        gradient_local[Nspatial_params + beta] += e_local*dbeta;
        gradient[Nspatial_params + beta] = dbeta;

    }

}

VMC* ASGD::minimize() {

    vmc->initialize();

    for (int cycle = 0; cycle < n_walkers * n_c; cycle++) {

        vmc->diffuse_walker();

        if ((cycle % n_c) == 0) {
            walkers[k] = vmc->clone_walker(vmc->original_walker);
            k++;
        }
    }


    for (int sample = 0; sample < SGDsamples; sample++) {

        E = 0;
        gradient = zeros(1, Nparams);
        gradient_local = zeros(1, Nparams);

        for (int k = 0; k < n_walkers; k++) {
            for (int cycle = 0; cycle < n_c_SGD; cycle++) {

                vmc->original_walker = vmc->clone_walker(walkers[k]);
                vmc->trial_walker = vmc->clone_walker(walkers[k]);

                vmc->diffuse_walker();

                vmc->calculate_energy_necessities(vmc->original_walker);
                double e_local = vmc->calculate_local_energy(vmc->original_walker);
                E += e_local;

                get_variational_gradients(e_local);

            }
        }


        int scale = n_walkers*n_c_SGD;

        E /= scale;
        gradient_tot = 2 * (gradient_local - gradient * E) / scale;

        double x = -dot(gradient_tot, gradient_old);

        t = (t_prev + f(x))*(t < 0);

        //output progress


        for (int param = 0; param < Nparams; param++) {

            step = a / (t + A) * gradient_tot[param];
            
            if (step * step > max_step * max_step) {
                step *= max_step / fabs(step);
            }
            
            if (param < Nspatial_params) {
                
                double alpha = vmc->get_orbitals_ptr()->get_parameter(param);
                vmc->get_orbitals_ptr()->set_parameter(abs(alpha - step), param);
                
            } else {
                
                double beta = vmc->get_jastrow_ptr()->get_parameter(param - Nspatial_params);
                vmc->get_jastrow_ptr()->set_parameter(abs(beta - step), param);
                
            }
            
            
            t_prev = t;
            gradient_old = gradient_tot;

        }
    }
    
    output();
    
    return vmc;
}