#include <sys/time.h>
#include <iostream>
#include <math.h>
#include <fstream>

#include "DMC.h"

using namespace std;

DMC::DMC() {

}

DMC::~DMC() {
    int i;

    delete system;

    for (i = 0; i < K * n_w_orig; i++) {
        delete Angry_mob[i];
    }

    delete[] Angry_mob;

}

DMC_BF::~DMC_BF() {
    int i;

    delete system;

    for (i = 0; i < K * n_w_orig; i++) {
        delete Angry_mob[i];
    }

    delete[] Angry_mob;

    free_matrix((void**) tmp_r);
}

int DMC::get_N_P() const {
    return n_p;
}

int DMC::get_dim() const {
    return dim;
}

bool DMC::singular(int k) {
    int i, j;
    double eps;

    eps = 0.1;

    for (i = 0; i < n_p; i++) {
        for (j = i + 1; j < n_p; j++) {
            if (Angry_mob[k]->r_rel[i][j] < eps) {
                return true;
            }
        }
    }

    return false;
}

void DMC::initialize_walkers() {
    int i, j, k;
    ifstream dist;
    dist.open("dist_out.dat");
    double pos;
    //    for (i=0; i < 100; i++){
    //        dist >> pos;
    //        cout << pos << endl;
    //    }
    //    exit(1);

    for (k = 0; k < n_w; k++) {
        Angry_mob[k] = new Walker(this);


        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                //Angry_mob[k]->r[i][j] = get_new_pos();
                dist >> pos;
                Angry_mob[k]->r[i][j] = pos;
            }
        }

        Angry_mob[k]->make_rel_matrix();


        //        while (singular(k)) {
        //
        //            for (i = 0; i < n_p; i++) {
        //                for (j = 0; j < dim; j++) {
        //                    Angry_mob[k]->r[i][j] = get_new_pos();
        //                }
        //            }
        //
        //            Angry_mob[k]->make_rel_matrix();
        //
        //        }

        if (is_importance()) {
            Angry_mob[k]->get_qforce(system);
        }

    }

    for (k = n_w; k < K * n_w_orig; k++) {
        Angry_mob[k] = new Walker(this);
    }
    dist.close();
}

void DMC::calc_gs_statistics() {
    //    int i, n;
    //    double N;
    //    bool tmp = true;
    //    n = n_c / 2;
    //
    //    N = 1;
    //    E = E_T;
    //    for (i = 1; i < n_c; i++) {
    //        if (population[i] == 0) {// || population[i + n] == 0) {
    //            if (tmp) {
    //                cout << "walkers died out." << endl;
    //                tmp = false;
    //            }
    //        }
    //        //        if (population[i] == 0 || population[i + n] == 0) {
    //        //
    //        //        } else {
    //        //            //N = N * (double) population[i] / population[i + n];
    //        //            
    //        //            E += 0.001 * log((double) population[0] / population[i]);
    //        //        }
    //    }
    //
    //    //E = E_T + 1 / (n * n * dt) * log(N);
    //    E /= (n_c - 1);


    cout << E / n_c << endl;
    //    cout << test / test2 << endl;
}

void DMC::Evolve_walker(int k) {
    int n;
    int branch_mean;

    branch_mean = int(GB + ran2(&ranseed)); //random int with mean=GB

    if (branch_mean > 1) {
        for (n = 0; n < branch_mean - 1; n++) {

            n_w += 1;
            copy_walker(k, n_w - 1);
        }

    } else if (branch_mean == 0) {
        Angry_mob[k]->is_murdered = true;
    }
}

void DMC::bury_the_dead() {
    int i, j, newborn, index, index2;
    int k = 0;
    //    for (i=0; i < n_w_now; i++){
    //        cout << Angry_mob[i]->is_murdered;
    //        if (Angry_mob[i]->is_murdered){
    //            cout << " dead ";
    //        } else {
    //            cout << " alive ";
    //        }
    //        k += 1;
    //        if (k == 10){
    //            cout << endl;
    //            k = 0;
    //        }
    //    }

    //    //TEST
    //    n_w = 6;
    //    n_w_now = 4;
    //    for (i=0; i < n_w; i++){
    //        delete Angry_mob[i];
    //        Angry_mob[i] = new Walker(this);
    //    }
    //    
    //    Angry_mob[0]->is_murdered = Angry_mob[1]->is_murdered = Angry_mob[3]->is_murdered = Angry_mob[4]->is_murdered = Angry_mob[5]->is_murdered = false;
    //    Angry_mob[2]->is_murdered = true;
    //    //



    newborn = n_w - n_w_now;
    j = 0;

    //    cout << "[";
    //    for (i = 0; i < n_w; i++) {
    //        if (Angry_mob[i]->is_murdered) {
    //            cout << "0, ";
    //        } else {
    //            cout << "1, ";
    //        }
    //    }
    //    cout << "]" << endl;

    for (i = 0; i < n_w; i++) {

        if (Angry_mob[i]->is_murdered) {


            //Fills the newborn into dead peoples homes
            if (j < newborn) {
                index = newborn + n_w_now - 1 - j;

                copy_walker(index, i);
                delete Angry_mob[index];
                Angry_mob[index] = new Walker(this);

                j++;
            } else {

                index2 = n_w_now - 1;

                //Find the last alive walker
                while (Angry_mob[index2]->is_murdered) {
                    index2--;
                    if (index2 < 0) {
                        cout << "no alive walkers after " << step << " cycles." << endl;
                    }
                }


                if (index2 > i) {
                    copy_walker(index2, i);
                    Angry_mob[index2]->is_murdered = true;
                }

            }




        }
    }

    //Finding the last alive walker below n_w_now
    for (i = 0; i < n_w_now; i++) {
        if (Angry_mob[i]->is_murdered) {

            //clearing up if someone is dead;
            delete Angry_mob[i];
            Angry_mob[i] = new Walker(this);
        } else {
            n_w = i + 1;
        }

    }

    //if there are more newborn than dead, at this point n_w = n_w_now;
    //adding the newborns which didn't get a home to that list:
    n_w += newborn - j;

    //    cout << "[";
    //    for (i = 0; i < n_w; i++) {
    //        if (Angry_mob[i]->is_murdered) {
    //            cout << "0, ";
    //        } else {
    //            cout << "1, ";
    //        }
    //    }
    //
    //    cout << "]" << endl;
    //    cout << n_w << endl;
    //
    //
}

void DMC::increase_walker_space() {
    int i, j, k, l;

    //creating tmp-replacement
    Walker **tmp = new Walker*[2 * K * n_w_orig];

    for (i = 0; i < 2 * K * n_w_orig; i++) {
        tmp[i] = new Walker(this);
    }

    //copying current walkers into the tmp
    for (l = 0; l < K * n_w_orig; l++) {
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                tmp[l]->r[i][j] = Angry_mob[l]->r[i][j];
            }

            for (k = i + 1; k < n_p; k++) {
                tmp[l]->r_rel[i][k] = tmp[l]->r_rel[k][i] = Angry_mob[l]->r[i][k];
            }
        }

        if (is_importance()) {
            for (i = 0; i < n_p; i++) {
                for (j = 0; j < dim; j++) {
                    tmp[l]->qforce[i][j] = Angry_mob[l]->qforce[i][j];
                }
            }
        }
        tmp[l]->is_murdered = Angry_mob[l]->is_murdered;
    }

    //deleting the old
    for (i = 0; i < K * n_w_orig; i++) {
        delete Angry_mob[i];
    }

    delete[] Angry_mob;

    //creating a new with more room
    K *= 2;
    Walker** Angry_mob = new Walker*[K * n_w_orig];

    for (i = 0; i < K * n_w_orig; i++) {
        Angry_mob[i] = new Walker(this);
    }



    //copying from tmp into the new
    for (l = 0; l < K * n_w_orig; l++) {
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                Angry_mob[l]->r[i][j] = tmp[l]->r[i][j];
            }

            for (k = i + 1; k < n_p; k++) {
                Angry_mob[l]->r_rel[i][k] = Angry_mob[l]->r_rel[k][i] = tmp[l]->r[i][k];
            }
        }

        if (is_importance()) {
            for (i = 0; i < n_p; i++) {
                for (j = 0; j < dim; j++) {
                    Angry_mob[l]->qforce[i][j] = tmp[l]->qforce[i][j];
                }
            }
        }
        Angry_mob[l]->is_murdered = tmp[l]->is_murdered;

    }


    //deleting tmp
    for (i = 0; i < K * n_w_orig; i++) {
        delete tmp[i];
    }

    delete[] tmp;

    this->Angry_mob = Angry_mob;

}

void DMC::copy_walker(int parent, int child) {
    int i, j, k;

    if (child >= K * n_w_orig) {
        cout << "NOT ENOUGH ROOM FOR WALKERS CYCLE: " << step << endl;
        cout << "INCREASING VALUE TO " << 2 * K * n_w_orig << endl;
        increase_walker_space();
    }

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            Angry_mob[child]->r[i][j] = Angry_mob[parent]->r[i][j];
        }

        for (k = i + 1; k < n_p; k++) {
            Angry_mob[child]->r_rel[i][k] = Angry_mob[child]->r_rel[k][i] = Angry_mob[parent]->r_rel[i][k];
        }
    }

    if (is_importance()) {
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                Angry_mob[child]->qforce[i][j] = Angry_mob[parent]->qforce[i][j];
            }
        }
    }

    Angry_mob[child]->is_murdered = Angry_mob[parent]->is_murdered;

}

void DMC_BF::update_pos(int k) {
    int i, j;

    

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            tmp_r[i][j] = Angry_mob[k]->r[i][j];
            Angry_mob[k]->r[i][j] += get_new_pos();
        }
    }

    Angry_mob[k]->make_rel_matrix();
    while (singular(k)) {
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                Angry_mob[k]->r[i][j] = tmp_r[i][j] + get_new_pos();
            }
        }

        Angry_mob[k]->make_rel_matrix();
    }

}

//////////////
//BRUTE FORCE
/////////////

DMC_BF::DMC_BF(double d, long RANSEED, int N_C, int N_W, System* SYSTEM, double DT, double e_t) {
    n_p = SYSTEM->get_N_P();
    dim = SYSTEM->get_dim();
    n_c = N_C;
    n_w = N_W;
    n_w_orig = N_W;
    D = d;
    ranseed = RANSEED;
    E_T = e_t;

    K = 2;

    system = SYSTEM;
    Angry_mob = new Walker*[K * N_W];

    dt = DT;
    std = sqrt(2 * DT * D);

    population = new int[n_c];
    tmp_r = (double **) matrix(n_p, dim, sizeof (double));
}

bool DMC_BF::is_importance() const {
    return false;
}

void DMC_BF::do_DMC() {
    int k, thresh;
    double V_x, V_y;

    ofstream file;
    file.open("e_vs_t.dat");

    initialize_walkers();
    E = E_T;
    thresh = n_c/1000;


    for (step = 0; step < n_c; step++) {

        n_w_now = n_w;
        population[step] = n_w;

        for (k = 0; k < n_w_now; k++) {

            V_x = system->get_pot_E(Angry_mob[k]);

            update_pos(k);

            V_y = system->get_pot_E(Angry_mob[k]);

            //Calculating the brancing G.F.
            GB = exp(-(0.5 * (V_y + V_x) - E_T) * dt);


            Evolve_walker(k);

        }

        bury_the_dead();



        E_T = E / (step + 1) + (1. / 2) * log((double) population[0] / n_w);
        E += E_T;

        if (step % thresh == 0) {
            cout << "\r";
            printf("%1.5f %1.5f %1.5f %1.5f%%",E_T, E / (step + 1), (double) n_w / population[0], (double) (step+1)/n_c*100);
            //cout << "\r" <<E_T << "\t" << E / (step + 1) << "\t" << (double) n_w / population[0] << "\t"<< (double) (step+1)/n_c*100<<"%"<< flush;
            file << (step+1)*dt <<"\t"<< E / (step + 1) << "\t" << (double) n_w / population[0] << endl;
        }
    }

    cout << "\r";
    cout << "Last Energy: " <<E / n_c << endl;
    double E_x = 2;
    double E_T0 = 2.00575;
    cout << "Improvement ratio (old error / new error): "<<abs((E_T0-E_x)/(E/n_c-E_x)) << endl;
    //calc_gs_statistics();
    file.close();

}

double DMC_BF::get_new_pos() {
    return gaussian_deviate(&ranseed) * std;
};









//////////////
//IMPORTANCE SAMPLING 
/////////////

DMC_Importance::DMC_Importance() {

}


