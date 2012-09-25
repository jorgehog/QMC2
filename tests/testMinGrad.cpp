#include "../src/QMCheaders.h"


void diffuse_walker(Walker* walkerPre, Walker* walkerNew){
    
}

void DoMin(){
    n_walkers = 10;
    Walker** walkers = new Walker*[n_walkers];
    Walker** trial_walkers = new Walker*[n_walkers];
    thermalization = 1E5;
    

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
                    DEBAG << E / ((cycle + 1)*(k + 1)) << endl;
                    get_variational_gradients(walkers[k], e_local);
                }

            }
            DEBAG << "Changed Walker" << endl;
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
            if (fabs(step) > max_step) {
                step *= max_step / fabs(step);
            }
            file << fabs(step) << "\t";
            double alpha = vmc->get_orbitals_ptr()->get_parameter(param);
            file << fabs(alpha - step) << "\t";
            vmc->get_orbitals_ptr()->set_parameter(fabs(alpha - step), param);
        }

        for (int param = 0; param < Njastrow_params; param++) {

            step = a / (t + A) * gradient_tot[Nspatial_params + param];
            if (step * step > max_step * max_step) {
                step *= max_step / fabs(step);
            }

            double beta = vmc->get_jastrow_ptr()->get_parameter(param);
            file << fabs(beta - step) << "\t";
            vmc->get_jastrow_ptr()->set_parameter(fabs(beta - step), param);


        }

        file << E << "\t";
        file << sumE / debug1 << endl;
        debug1++;

        t_prev = t;
        gradient_old = gradient_tot;

    }

    cout << "done" << endl;
}


int main(){
    
}