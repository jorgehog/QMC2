/* 
 * File:   ASGD.cpp
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#include "../../QMCheaders.h"

using namespace arma;

ASGD::ASGD(VMC* vmc, MinimizerParams & mP)
: Minimizer(vmc, mP.alpha, mP.beta) {

    this->n_c = mP.n_cm;
    this->thermalization = mP.thermalization;
    this->n_c_SGD = mP.n_c_SGD;
    this->SGDsamples = mP.SGDsamples;
    this->n_walkers = mP.n_walkers;

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

    t_prev = A;
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
//    ofstream file;
//    file.open("alpha.dat");
//    DEBAG.open("sammenliknnBF_Num.dat");
    int debug1;
    double aGrad, bGrad, sumE;
    aGrad = bGrad = sumE = 0;
    debug1 = 1;
    //


    vmc->initialize();

    int k = 0;
    for (int cycle = 1; cycle <= thermalization + n_walkers * n_c; cycle++) {
        vmc->diffuse_walker(vmc->original_walker, vmc->trial_walker);

        if ((cycle > thermalization) && (cycle % n_c) == 0) {
            vmc->copy_walker(vmc->original_walker, walkers[k]);
            vmc->copy_walker(vmc->trial_walker, trial_walkers[k]);
            k++;
        }
    }

    //D:100% match this far for BF and IS. Walker copy funker for alle.

    int corrLength = 1;
    for (int sample = 1; sample <= SGDsamples; sample++) {

        E = 0;
        gradient = zeros(1, Nparams);
        gradient_local = zeros(1, Nparams);

        for (k = 0; k < n_walkers; k++) {
            for (int cycle = 0; cycle < n_c_SGD; cycle++) {

                vmc->diffuse_walker(walkers[k], trial_walkers[k]);

                if (cycle % corrLength == 0) {
                    vmc->calculate_energy_necessities(walkers[k]);
                    double e_local = vmc->calculate_local_energy(walkers[k]);
                    E += e_local;
//                    DEBAG << E / ((cycle + 1)*(k + 1)) << endl;
                    get_variational_gradients(walkers[k], e_local);
                }

            }
//            DEBAG << "Changed Walker" << endl;
        }

        int scale = n_walkers * n_c_SGD / corrLength;

        //        cout << E / scale << endl;
        E /= scale;

        gradient_tot = 2 * (gradient_local - gradient * E) / scale;

        double x = -dot(gradient_tot, gradient_old);

        t = t_prev + f(x);
        if (t < 0) {
            t = 0;
        }

        //        output progress
        if ((sample % 1) == 0) {
            output("cycle:", sample);
            cout << gradient_tot << endl;
            cout << E << endl;
        }


        //        
        aGrad += gradient_tot[0];
        bGrad += gradient_tot[1];
        sumE += E;
//        file << aGrad / debug1 << "\t" << bGrad / debug1 << "\t";
        //


        for (int param = 0; param < Nspatial_params; param++) {

            step = a / (t + A) * gradient_tot[param];
            if (fabs(step) > max_step) {
                step *= max_step / fabs(step);
            }
//            file << fabs(step) << "\t";
            double alpha = vmc->get_orbitals_ptr()->get_parameter(param);
//            file << fabs(alpha - step) << "\t";
            vmc->get_orbitals_ptr()->set_parameter(fabs(alpha - step), param);
        }

        for (int param = 0; param < Njastrow_params; param++) {

            step = a / (t + A) * gradient_tot[Nspatial_params + param];
            if (step * step > max_step * max_step) {
                step *= max_step / fabs(step);
            }

            double beta = vmc->get_jastrow_ptr()->get_parameter(param);
//            file << fabs(beta - step) << "\t";
            vmc->get_jastrow_ptr()->set_parameter(fabs(beta - step), param);


        }

//        file << E << "\t";
//        file << sumE / debug1 << endl;
        debug1++;

        t_prev = t;
        gradient_old = gradient_tot;

    }

    output("Finished minimizing. Final parameters:", -1);
//    file.close();
//    DEBAG.close();
    vmc->accepted = 0;
    return vmc;
}

double ASGD::TESTWF(Walker* walker) {
    double alpha = vmc->get_orbitals_ptr()->get_parameter(0);
    double X = walker->r(0, 0);
    return exp(-alpha * alpha * X * X);
}

double ASGD::TEST_E(Walker* walker) {
    double alpha = vmc->get_orbitals_ptr()->get_parameter(0);
    double X = walker->r(0, 0);
    return alpha * alpha + X * X * (0.5 - 2 * alpha * alpha*alpha*alpha);
}

double ASGD::TEST_G(Walker* walker_post, Walker* walker_pre) {
    double g_ratio = 0;
    double D = 0.5;
    double timestep = 0.05;
    int particle = 0;
    int j = 0;

    g_ratio += 0.5 * (walker_pre->qforce(particle, j) + walker_post->qforce(particle, j))*
            (D * timestep * 0.5 * (walker_pre->qforce(particle, j) - walker_post->qforce(particle, j))
            - walker_post->r(particle, j) + walker_pre->r(particle, j));


    return exp(g_ratio);
}

void ASGD::TEST_DIFF(Walker* original, Walker* trial) {
    double ALP = vmc->get_orbitals_ptr()->get_parameter(0);
    double X = original->r(0, 0);
    double dt = 0.05;
    double std = sqrt(dt);

    double QF = -2 * ALP * ALP*X;
    original->qforce(0, 0) = QF;

    trial->r(0, 0) = X + 0.5 * dt * QF + gaussian_deviate(&random_seed) * std;
    double Xnew = trial->r(0, 0);
    trial->qforce(0, 0) = -2 * ALP * ALP*Xnew;

    double A = TESTWF(trial) / TESTWF(original);
    A = A * A * TEST_G(trial, original);

    if (ran3(&random_seed) <= A) {
        original->r(0, 0) = trial->r(0, 0);
    } 
}

VMC* ASGD::minimizeTEST() {

    int dim = 1;
    int n_p = 1;




    random_seed = -1.234;
    double dt = 0.05;





    double std = sqrt(dt);

    vmc->original_walker->r(0, 0) = gaussian_deviate(&random_seed) * std;
    vmc->trial_walker->r(0, 0) = vmc->original_walker->r(0, 0);

    int k = 0;
    for (int cycle = 1; cycle <= thermalization + n_walkers * n_c; cycle++) {

        TEST_DIFF(vmc->original_walker, vmc->trial_walker);

        if ((cycle > thermalization) && (cycle % n_c) == 0) {
            vmc->copy_walker(vmc->original_walker, walkers[k]);
            vmc->copy_walker(vmc->trial_walker, trial_walkers[k]);
            k++;
        }
    }

    //D:100% match this far for BF and IS. Walker copy funker for alle.
    double gradient_old = 0;
    double gradient_tot = 0;
    double ALPH;
    int corrLength = 1;
    for (int sample = 1; sample <= SGDsamples; sample++) {

        E = 0;
        double gradient = 0;
        double gradient_local = 0;

        for (k = 0; k < n_walkers; k++) {
            for (int cycle = 0; cycle < n_c_SGD; cycle++) {

                TEST_DIFF(walkers[k], trial_walkers[k]);



                double e_local = TEST_E(walkers[k]);
                E += e_local;

                ALPH = vmc->get_orbitals_ptr()->get_parameter(0);
                double X = walkers[k]->r(0, 0);
                double dalpha = -ALPH * X*X;

                gradient_local += e_local * dalpha;
                gradient += dalpha;


            }

        }

        int scale = n_walkers * n_c_SGD / corrLength;

        //        cout << E / scale << endl;
        E /= scale;
//        cout << gradient_local << " " << gradient*E << endl;
        gradient_tot = 2 * (gradient_local - gradient * E) / scale;

        double x = -gradient_tot*gradient_old;

        t = t_prev + f(x);
        if (t < 0) {
            t = 0;
        }

        //        output progress
        if ((sample % 100) == 0) {
            cout << "CYCLE " << sample << endl;
            ALPH = vmc->get_orbitals_ptr()->get_parameter(0);
            cout << ALPH << endl;
            cout << E << " CF: " << 0.5 * ALPH * ALPH + 1.0 / (8 * ALPH * ALPH) << endl;
            cout << gradient_tot << " CF: " << ALPH - 1.0 / (4 * ALPH * ALPH * ALPH) << endl;
        }






        step = a / (t + A) * gradient_tot;
        if (fabs(step) > max_step) {
            step *= max_step / fabs(step);
        }

        double alpha = vmc->get_orbitals_ptr()->get_parameter(0);

        vmc->get_orbitals_ptr()->set_parameter(fabs(alpha - step), 0);



        t_prev = t;
        gradient_old = gradient_tot;

    }
    ALPH = vmc->get_orbitals_ptr()->get_parameter(0);
    cout << "FOUND" << ALPH << " CF: " << 1. / sqrt(2) << endl;

    return vmc;
}