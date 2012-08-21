#include "vmc_old.h"
#include "lib.h"
#include <sys/time.h>
#include <iostream>
#include <math.h>
#include <fstream>



using namespace std;

/*
 CONSTRUCTORS
 */


Pade_Jastrow::Pade_Jastrow(const VMC& vmc, double BETA) {
    n_p = vmc.n_p;
    n2 = vmc.n2;
    dim = vmc.dim;
    beta = BETA;
    a = (double **) matrix(n_p, n_p, sizeof (double));

}

void Pade_Jastrow::clear() {
    free_matrix((void **) a);
}

qdots::qdots(const VMC &vmc, const Harmonic_osc &pot, double ALPHA) {
    n_p = vmc.n_p;
    n2 = vmc.n2;
    dim = vmc.dim;

    alpha = ALPHA;
    w = pot.w;

    a_sym = 1. / 3;
    a_asym = 1.0;

    inv_new = (double **) matrix(n2, n_p, sizeof (double));
    inv_old = (double **) matrix(n2, n_p, sizeof (double));

    s_down = (double **) matrix(n2, n2, sizeof (double));
    s_up = (double **) matrix(n2, n2, sizeof (double));
}

Wavefunction::Wavefunction(const VMC& vmc) {
    n_p = vmc.n_p;
    n2 = vmc.n2;
    dim = vmc.dim;



    r = (double **) matrix(n_p, dim, sizeof (double));
    r_rel = (double **) matrix(n_p, n_p, sizeof (double));

    qforce = (double **) matrix(n_p, dim, sizeof (double));

    spatial_grad = (double **) matrix(n_p, dim, sizeof (double));
    jast_grad = (double **) matrix(n_p, dim, sizeof (double));

}

Wavefunction::Wavefunction() {

}

Kinetics_cf::Kinetics_cf(VMC &vmc) {
    n_p = vmc.n_p;
    dim = vmc.dim;

    cf = true;
}

Kinetics_num::Kinetics_num(VMC &vmc) {
    n_p = vmc.n_p;
    dim = vmc.dim;

    h = 0.001;
    h2 = 1000000;

    wfplus = Wavefunction(vmc);
    wfminus = Wavefunction(vmc);

    cf = false;
}

Kinetics_num::Kinetics_num(VMC &vmc, double H) {
    n_p = vmc.n_p;
    dim = vmc.dim;

    h = H;
    h2 = 1 / (H * H);

    wfplus = Wavefunction(vmc);
    wfminus = Wavefunction(vmc);
}

Harmonic_osc::Harmonic_osc(const VMC& vmc, double W) {
    n_p = vmc.n_p;
    dim = vmc.dim;
    w = W;
}

Importance::Importance(int N_C, int N_P, int DIM, double TIMESTEP) {
    n_c = N_C;
    n_p = N_P;
    n2 = n_p / 2;
    dim = DIM;

    timestep = TIMESTEP;
    sqrt_timestep = sqrt(timestep);

    local_e = local_e2 = 0;

    random_seed = -time(NULL);
    map_wf = false;

    D = 0.5;
}

Brute_Force::Brute_Force(int N_C, int N_P, int DIM, double STEPLENGTH) {
    n_c = N_C;
    n_p = N_P;
    n2 = n_p / 2;
    dim = DIM;
    steplength = STEPLENGTH;

    local_e = local_e2 = 0;
    random_seed = -time(NULL);
}

/*
 MEMBER FUNCTIONS
 */

double Kinetics_cf::get_KE(Wavefunction& wf, System& system, Jastrow_factor &jastrow) {
    int i, j;
    double xterm, e_kinetic;


    e_kinetic = xterm = 0;

    //the X-term
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            xterm += wf.jast_grad[i][j] * wf.spatial_grad[i][j];
        }
    }
    e_kinetic += 2 * xterm;

    //laplacian terms
    e_kinetic += system.get_lapl_sum(wf) + jastrow.get_lapl_sum(wf);


    return -0.5 * e_kinetic;
}

double Kinetics_num::get_KE(Wavefunction& wf, System& system, Jastrow_factor &jastrow) {
    int i, j;
    double e_kinetic;

    //kinetic energy:
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wfplus.r[i][j] = wfminus.r[i][j] = wf.r[i][j];
        }
    }

    e_kinetic = 0;
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wfplus.r[i][j] = wf.r[i][j] + h;
            wfminus.r[i][j] = wf.r[i][j] - h;

            wfplus.make_rel_matrix();
            wfminus.make_rel_matrix();

            system.get_wf_val(wfplus, jastrow);
            system.get_wf_val(wfminus, jastrow);

            e_kinetic -= (wfminus.wf_val + wfplus.wf_val - 2 * wf.wf_val);
            wfplus.r[i][j] = wf.r[i][j];
            wfminus.r[i][j] = wf.r[i][j];
        }
    }

    e_kinetic = 0.5 * h2 * e_kinetic / wf.wf_val;

    return e_kinetic;
}

void Kinetics_cf::get_QF(Wavefunction& wf, System &system, Jastrow_factor &jastrow) {
    int i, j;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wf.qforce[i][j] = 2 * (wf.jast_grad[i][j] + wf.spatial_grad[i][j]);
        }
    }
}

void Kinetics_num::get_QF(Wavefunction& wf, System &system, Jastrow_factor &jastrow) {
    int i, j;


    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wfplus.r[i][j] = wfminus.r[i][j] = wf.r[i][j];
        }
    }


    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wfplus.r[i][j] = wf.r[i][j] + h;
            wfminus.r[i][j] = wf.r[i][j] - h;

            wfplus.make_rel_matrix();
            wfminus.make_rel_matrix();

            system.get_wf_val(wfplus, jastrow);
            system.get_wf_val(wfminus, jastrow);

            wf.qforce[i][j] = (wfplus.wf_val - wfminus.wf_val) / (wf.wf_val * h);

            wfplus.r[i][j] = wf.r[i][j];
            wfminus.r[i][j] = wf.r[i][j];

        }
    }
}

void Kinetics_cf::clear() {

}

void Kinetics_num::clear() {
    wfminus.clear();
    wfplus.clear();
}

double Harmonic_osc::get_pot_E(const Wavefunction & wf) const {
    int i;
    double e_potential;

    e_potential = 0;

    // contribution from oscillator part 
    for (i = 0; i < n_p; i++) {
        e_potential += 0.5 * w * w * wf.get_r_i2(i);
    }

    return e_potential;
}

double Potential::get_Coulomb(const Wavefunction& wf) const {
    int i, j;
    double e_coulomb;

    e_coulomb = 0;
    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            e_coulomb += 1 / wf.r_rel[i][j];
        }
    }

    return e_coulomb;
}

void Jastrow_factor::get_grad(Wavefunction & wf) {
    int i, k;

    for (i = 0; i < n_p; i++) {
        for (k = 0; k < dim; k++) {
            wf.jast_grad[i][k] = get_deriv(wf, i, k);
        }
    }
}

double Pade_Jastrow::get_lapl_sum(const Wavefunction & wf) const {
    int k, j, d;
    double sum1, sum2, b_kj;

    sum1 = sum2 = 0;

    for (k = 0; k < n_p; k++) {
        for (j = k + 1; j < n_p; j++) {
            b_kj = 1 + beta * wf.r_rel[k][j];
            sum2 += a[k][j] * (1 - beta * wf.r_rel[k][j]) / (wf.r_rel[k][j] * b_kj * b_kj * b_kj);
        }

        for (d = 0; d < dim; d++) {
            sum1 += wf.jast_grad[k][d] * wf.jast_grad[k][d];
        }
    }

    return sum1 + 2 * sum2;
}

double Pade_Jastrow::get_deriv(Wavefunction& wf, int i, int k) {
    int j;
    double b_ij, deriv;

    deriv = 0;
    for (j = 0; j < i; j++) {
        b_ij = 1.0 + beta * wf.r_rel[i][j];
        deriv += a[i][j] * (wf.r[i][k] - wf.r[j][k]) / (wf.r_rel[i][j] * b_ij * b_ij);
    }
    for (j = i + 1; j < n_p; j++) {
        b_ij = 1.0 + beta * wf.r_rel[i][j];
        deriv += a[i][j] * (wf.r[i][k] - wf.r[j][k]) / (wf.r_rel[i][j] * b_ij * b_ij);
    }


    return deriv;
}

double Pade_Jastrow::get_j_ratio(Wavefunction &wf_new, Wavefunction &wf_old, int i) {
    int j;
    double j_ratio;

    j_ratio = 0;
    for (j = 0; j < n_p; j++) {
        j_ratio += a[i][j] * (wf_new.r_rel[i][j] / (1.0 + beta * wf_new.r_rel[i][j]) -
                wf_old.r_rel[i][j] / (1.0 + beta * wf_old.r_rel[i][j]));
    }


    return exp(j_ratio);
}

void Pade_Jastrow::initialize_params(const System & system) const {
    int i, j;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < n_p; j++) {

            if ((j < n2 && i < n2) || (j >= n2 && i >= n2)) {
                a[i][j] = system.a_sym;
            } else {
                a[i][j] = system.a_asym;
            }
        }
    }
}

void fermions::invert_slaters() {
    int i, j;
    double d;
    // allocate space in memory
    int* indx = new int[n2];
    double* col = new double[n2];
    double** y = (double **) matrix(n2, n2, sizeof (double));


    //SPIN UP:
    ludcmp(s_up, n2, indx, &d); // LU decompose a

    // find inverse of a by columns
    for (j = 0; j < n2; j++) {
        // initialize right-side of linear equations
        for (i = 0; i < n2; i++) col[i] = 0.0;
        col[j] = 1.0;
        lubksb(s_up, n2, indx, col);
        // save result in y
        for (i = 0; i < n2; i++) y[i][j] = col[i];
    } //j-loop over columns
    // return the inverse matrix in a
    for (i = 0; i < n2; i++) {
        for (j = 0; j < n2; j++) s_up[i][j] = y[i][j];
    }


    //SPIN DOWN:
    ludcmp(s_down, n2, indx, &d); // LU decompose a

    // find inverse of a by columns
    for (j = 0; j < n2; j++) {
        // initialize right-side of linear equations
        for (i = 0; i < n2; i++) col[i] = 0.0;
        col[j] = 1.0;
        lubksb(s_down, n2, indx, col);
        // save result in y
        for (i = 0; i < n2; i++) y[i][j] = col[i];
    } //j-loop over columns
    // return the inverse matrix in a
    for (i = 0; i < n2; i++) {
        for (j = 0; j < n2; j++) s_down[i][j] = y[i][j];
    }

    delete[] col;
    delete[] indx;
    free_matrix((void **) y);

}

double fermions::get_spatial_ratio(const Wavefunction& wf_new, int i) {
    int q_num;

    s_ratio = 0;
    for (q_num = 0; q_num < n2; q_num++) {
        s_ratio += phi(wf_new, i, q_num) * inv_old[q_num][i];
    }

    return s_ratio;
}

void fermions::initialize(Wavefunction& wf, Jastrow_factor &jastrow, Kinetics &kinetics) {
    make_merged_inv(wf);
    if (kinetics.cf) {
        get_init_grad(wf, jastrow);
    } else {
        get_wf_val(wf, jastrow);
    }
}

void fermions::initialize_slaters(const Wavefunction & wf) {
    int i, q_num;

    for (i = 0; i < n2; i++) {
        for (q_num = 0; q_num < n2; q_num++) {
            s_up[i][q_num] = phi(wf, i, q_num);
            s_down[i][q_num] = phi(wf, i + n2, q_num);
        }
    }
}

void fermions::make_merged_inv(const Wavefunction & wf) {
    int i, j;

    initialize_slaters(wf);
    // initializing the spin slaters  
    //    for (i = 0; i < n2; i++) {
    //        for (q_num = 0; q_num < n2; q_num++) {
    //            s_up[i][q_num] = phi(wf, i, q_num);
    //            s_down[i][q_num] = phi(wf, i + n2, q_num);
    //        }
    //    }

    invert_slaters();

    //merging the inverse matrices
    for (i = 0; i < n2; i++) {
        for (j = 0; j < n2; j++) {
            inv_old[i][j] = s_up[i][j];
            inv_old[i][j + n2] = s_down[i][j];
        }
    }



}

void fermions::get_init_grad(Wavefunction& wf, Jastrow_factor & jastrow) {
    int i, j, k;

    for (i = 0; i < n_p; i++) {
        for (k = 0; k < dim; k++) {
            wf.spatial_grad[i][k] = 0;
            wf.jast_grad[i][k] = jastrow.get_deriv(wf, i, k);

            for (j = 0; j < n2; j++) {
                wf.spatial_grad[i][k] += del_phi(wf, i, j, k) * inv_old[j][i];
            }
        }
    }
}

void fermions::get_new_grad(const Wavefunction& wf_old, Wavefunction &wf_new, int particle) {
    int i, j, k, start;


    if (particle >= n2) {
        start = n2;
    } else {
        start = 0;
    }

    //if the last metropolis-test went through, there is no need to overwrite data
    if (accepted_last == false) {
        for (i = n2 - start; i < n_p - start; i++) {
            for (k = 0; k < dim; k++) {
                wf_new.spatial_grad[i][k] = wf_old.spatial_grad[i][k];
            }
        }
    }


    //updating the part with the same spin as the moved particle
    for (i = start; i < n2 + start; i++) {
        for (k = 0; k < dim; k++) {
            wf_new.spatial_grad[i][k] = 0;

            for (j = 0; j < n2; j++) {
                wf_new.spatial_grad[i][k] += del_phi(wf_new, i, j, k) * inv_new[j][i];
            }

        }
    }
}

void fermions::calc_for_newpos(const Wavefunction& wf_old, Wavefunction& wf_new, Kinetics &kinetics, Jastrow_factor &jastrow, int i) {
    int j, k, l, start;
    double sum;

    //UPDATING THE INVERSE

    if (i >= n2) {
        start = n2;
    } else {
        start = 0;
    }

    //if the last metropolis-test went through, there is no need to overwrite data
    if (accepted_last == false) {
        for (j = n2 - start; j < n_p - start; j++) {
            for (k = 0; k < n2; k++) {
                inv_new[k][j] = inv_old[k][j];
            }
        }
    }

    //updating the part of the inverse with the same spin as particle i
    for (k = 0; k < n2; k++) {
        for (j = start; j < n2 + start; j++) {
            sum = 0;
            if (j == i) {
                for (l = 0; l < n2; l++) {
                    sum += phi(wf_old, i, l) * inv_old[l][j];
                }
                inv_new[k][j] = inv_old[k][i] / s_ratio*sum;
            } else {
                for (l = 0; l < n2; l++) {
                    sum += phi(wf_new, i, l) * inv_old[l][j];
                }
                inv_new[k][j] = inv_old[k][j] - inv_old[k][i] / s_ratio*sum;
            }
        }
    }

    //GETTING NEW VALS FOR QFORCE
    if (kinetics.cf) {
        get_new_grad(wf_old, wf_new, i);
        jastrow.get_grad(wf_new);
    } else {
        get_wf_val(wf_new, jastrow);
    }

}

void fermions::update_old(int i) {
    int j, k, start;

    if (i >= n2) {
        start = n2;
    } else {
        start = 0;
    }

    //updating the parts with the same spin as the moved particle
    for (k = start; k < start + n2; k++) {
        for (j = 0; j < n2; j++) {
            inv_old[j][k] = inv_new[j][k];
        }
    }
}

void fermions::clear() {
    free_matrix((void **) inv_new);
    free_matrix((void **) inv_old);
    free_matrix((void **) s_up);
    free_matrix((void **) s_down);
}

double qdots::phi(const Wavefunction& wf, int particle, int q_num) {
    double r2, H, x, y;

    r2 = wf.get_r_i2(particle);
    x = wf.r[particle][0];
    y = wf.r[particle][1];


    if (q_num == 0) {
        H = 1;
    } else if (q_num == 1) {
        H = 2 * x;
    } else if (q_num == 2) {
        H = 2 * y;
    } else if (q_num == 3) {
        H = 4 * x * x - 2;
    } else if (q_num == 4) {
        H = 4 * y * y - 2;
    } else if (q_num == 5) {
        H = 4 * x*y;
    } else if (q_num == 6) {
        H = 4 * x * (2 * y * y - 1);
    } else if (q_num == 7) {
        H = 4 * y * (2 * x * x - 1);
    } else if (q_num == 8) {
        H = 8 * y * y * y - 12 * y;
    } else if (q_num == 9) {
        H = 8 * x * x * x - 12 * x;
    } else if (q_num == 10) {
        H = 8 * y * x * (2 * x * x - 3);
    } else if (q_num == 11) {
        H = 8 * y * x * (2 * y * y - 3);
    } else if (q_num == 12) {
        H = 16 * y * y * (y * y - 3) + 12;
    } else if (q_num == 13) {
        H = 16 * x * x * (x * x - 3) + 12;
    } else if (q_num == 14) {
        H = 4 * (2 * x * x - 1)*(2 * y * y - 1);
    } else {
        cout << "Mismatching quantum number: " << q_num << endl;
    }

    return H * exp(-0.5 * alpha * w * r2);
}

double qdots::del_phi(const Wavefunction& wf, int particle, int q_num, int d) {
    double r2, H, x, y;

    r2 = wf.get_r_i2(particle);

    x = wf.r[particle][0];
    y = wf.r[particle][1];


    //Hermite polynomials (up to 4)
    if (q_num == 0) {
        H = -w * alpha * wf.r[particle][d];
    } else if (q_num == 1) {
        if (d == 0) {
            H = 2 * (1 - alpha * w * x * x);
        } else {
            H = -2 * alpha * w * x*y;
        }
    } else if (q_num == 2) {
        if (d == 0) {
            H = -2 * alpha * w * y*x;
        } else {
            H = 2 * (1 - alpha * w * y * y);
        }

    } else if (q_num == 3) {
        if (d == 0) {
            H = 2 * x * (4 + alpha * w - 2 * x * x * alpha * w);
        } else {
            H = (4 * x * x - 2)*(-alpha * w * y);
        }
    } else if (q_num == 4) {
        if (d == 0) {
            H = (4 * y * y - 2)*(-alpha * w * x);
        } else {
            H = 2 * y * (4 + alpha * w - 2 * y * y * alpha * w);
        }
    } else if (q_num == 5) {
        if (d == 0) {
            H = 4 * y - 4 * x * x * y * w*alpha;
        } else {
            H = 4 * x - 4 * x * y * y * w*alpha;
        }
    } else if (q_num == 6) {
        if (d == 0) {
            H = -4 * (-1 + alpha * w * x * x)*(-1 + 2 * y * y);
        } else {
            H = -4 * x * y * (alpha * w * (2 * y * y - 1) - 4);
        }
    } else if (q_num == 7) {
        if (d == 0) {
            H = -4 * x * y * (alpha * w * (2 * x * x - 1) - 4);
        } else {
            H = -4 * (-1 + alpha * w * y * y)*(-1 + 2 * x * x);
        }
    } else if (q_num == 8) {
        if (d == 0) {
            H = 4 * alpha * w * y * (3 - 2 * y * y) * x;
        } else {
            H = 4 * (y * y * (alpha * w * (3 - 2 * y * y) + 6) - 3);
        }
    } else if (q_num == 9) {
        if (d == 0) {
            H = 4 * (x * x * (alpha * w * (3 - 2 * x * x) + 6) - 3);
        } else {
            H = 4 * alpha * w * x * (3 - 2 * x * x) * y;
        }
    } else if (q_num == 10) {
        if (d == 0) {
            H = -8 * y * (alpha * w * x * x * (2 * x * x - 3) - 3 * (2 * x * x - 1));
        } else {
            H = -8 * x * (2 * x * x - 3)*(alpha * w * y * y - 1);
        }
    } else if (q_num == 11) {
        if (d == 0) {
            H = -8 * y * (2 * y * y - 3)*(alpha * w * x * x - 1);
        } else {
            H = -8 * x * (alpha * w * y * y * (2 * y * y - 3) - 3 * (2 * y * y - 1));
        }
    } else if (q_num == 12) {
        if (d == 0) {
            H = -4 * alpha * w * x * (4 * y * y * y * y - 12 * y * y + 3);
        } else {
            H = -4 * y * (3 * (alpha * w + 8) + 4 * alpha * w * y * y * y * y - 4 * y * y * (3 * alpha * w + 4));
        }
    } else if (q_num == 13) {
        if (d == 0) {
            H = -4 * x * (3 * (alpha * w + 8) + 4 * alpha * w * x * x * x * x - 4 * x * x * (3 * alpha * w + 4));
        } else {
            H = -4 * alpha * w * y * (4 * x * x * x * x - 12 * x * x + 3);
        }
    } else if (q_num == 14) {
        if (d == 0) {
            H = -4 * x * (2 * y * y - 1)*(alpha * w * (2 * x * x - 1) - 4);
        } else {
            H = -4 * y * (2 * x * x - 1)*(alpha * w * (2 * y * y - 1) - 4);
        }
    } else {
        cout << "Mismatching quantum number: " << q_num << endl;
    }

    return H * exp(-0.5 * alpha * w * r2);
}

double qdots::lapl_phi(const Wavefunction& wf, int particle, int q_num) const {
    double r2, H, x, y, aw2;

    r2 = wf.get_r_i2(particle);
    x = wf.r[particle][0];
    y = wf.r[particle][1];
    aw2 = alpha * alpha * w * w;


    //Hermite polynomials (up to 4)
    if (q_num == 0) {
        H = w * alpha * (w * alpha * r2 - dim);
    } else if (q_num == 1) {
        H = 2 * x * w * alpha * (w * alpha * r2 - 4);
    } else if (q_num == 2) {
        H = 2 * y * w * alpha * (w * alpha * r2 - 4);
    } else if (q_num == 3) {
        H = 2 * (4 + w * alpha * (2 - 12 * x * x + w * alpha * r2 * (2 * x * x - 1)));
    } else if (q_num == 4) {
        H = 2 * (4 + w * alpha * (2 - 12 * y * y + w * alpha * r2 * (2 * y * y - 1)));
    } else if (q_num == 5) {
        H = 4 * alpha * w * x * y * (alpha * w * r2 - 6);
    } else if (q_num == 6) {
        H = 4 * x * (aw2 * r2 * (2 * y * y - 1) + 4 * alpha * w * (1 - 4 * y * y) + 4);
    } else if (q_num == 7) {
        H = 4 * y * (aw2 * r2 * (2 * x * x - 1) + 4 * alpha * w * (1 - 4 * x * x) + 4);
    } else if (q_num == 8) {
        H = 4 * y * (aw2 * r2 * (2 * y * y - 3) + 4 * alpha * w * (3 - 4 * y * y) + 12);
    } else if (q_num == 9) {
        H = 4 * x * (aw2 * r2 * (2 * x * x - 3) + 4 * alpha * w * (3 - 4 * x * x) + 12);
    } else if (q_num == 10) {
        H = 8 * x * y * (aw2 * r2 * (2 * x * x - 3) + 2 * alpha * w * (9 - 10 * x * x) + 12);
    } else if (q_num == 11) {
        H = 8 * x * y * (aw2 * r2 * (2 * y * y - 3) + 2 * alpha * w * (9 - 10 * y * y) + 12);
    } else if (q_num == 12) {
        H = 4 * (aw2 * r2 * (3 + 4 * y * y * (y * y - 3)) + 2 * alpha * w * (36 * y * y - 3 - 20 * y * y * y * y) + 24 * (y * y - 1));
    } else if (q_num == 13) {
        H = 4 * (aw2 * r2 * (3 + 4 * x * x * (x * x - 3)) + 2 * alpha * w * (36 * x * x - 3 - 20 * x * x * x * x) + 24 * (x * x - 1));
    } else if (q_num == 14) {
        H = 4 * (aw2 * r2 * (4 * x * x * y * y + 1 + 2 * r2) - 2 * alpha * w * (1 + 6 * r2 + 20 * x * x * y * y) + 8 * (r2 + 1));
    } else {
        cout << "Mismatching quantum number: " << q_num << endl;
    }

    return H * exp(-0.5 * alpha * w * r2);
}

double fermions::get_lapl_sum(const Wavefunction & wf) const {
    int i, j;
    double sum;

    sum = 0;
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < n_p / 2; j++) {
            sum += lapl_phi(wf, i, j) * inv_old[j][i];
        }
    }

    return sum;
}

double Wavefunction::get_r_i2(int i) const {
    int j;
    double r2;

    r2 = 0;
    for (j = 0; j < dim; j++) {
        r2 += r[i][j] * r[i][j];
    }

    return r2;
}

double Wavefunction::abs_relative(int i, int j) {
    int k;
    double r_ij, tmp;

    r_ij = 0;
    for (k = 0; k < dim; k++) {
        tmp = (r[i][k] - r[j][k]);
        r_ij += tmp*tmp;
    }
    r_ij = sqrt(r_ij);

    return r_ij;
}

void Wavefunction::make_rel_matrix() {
    int i, j;


    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            r_rel[i][j] = r_rel[j][i] = abs_relative(i, j);
        }
    }
}

void Wavefunction::update_old(Wavefunction& wf_new, Kinetics &kinetics, int i) {
    int j, k, start;

    for (j = 0; j < dim; j++) {
        r[i][j] = wf_new.r[i][j];
        for (k = 0; k < n_p; k++) {
            qforce[k][j] = wf_new.qforce[k][j];
        }
    }

    for (k = 0; k < n_p; k++) {
        r_rel[k][i] = r_rel[i][k] = wf_new.r_rel[k][i];
    }

    if (kinetics.cf) {

        if (i >= n2) {
            start = n2;
        } else {
            start = 0;
        }

        //updating the parts with the same spin as the moved particle
        for (k = start; k < start + n2; k++) {
            for (j = 0; j < dim; j++) {
                spatial_grad[k][j] = wf_new.spatial_grad[k][j];
            }
        }
        for (k = 0; k < n_p; k++) {
            for (j = 0; j < dim; j++) {
                jast_grad[k][j] = wf_new.jast_grad[k][j];
            }
        }
    } else {
        wf_val = wf_new.wf_val;
    }
}

void Wavefunction::clear() {
    free_matrix((void **) r);
    free_matrix((void **) r_rel);
    free_matrix((void **) qforce);
    free_matrix((void **) spatial_grad);
    free_matrix((void **) jast_grad);
}

bool VMC::metropolis_test() {
    random_number = ran2(&random_seed);
    if (random_number <= full_ratio) {
        return true;
    } else {
        return false;
    }
}

void VMC::output() {
    cout << endl;
    cout << "energy" << "\t\t" << "variance" << endl;
    cout << local_e << "\t\t" << local_e2 - local_e * local_e << endl;
    cout << endl;
}

void VMC::free_memory(Wavefunction& wf_old, Wavefunction& wf_new, Jastrow_factor& jastrow, System& system, Kinetics & kinetics) {
    wf_old.clear();
    wf_new.clear();
    jastrow.clear();
    system.clear();
    kinetics.clear();
}

void Importance::update_walker(Wavefunction& wf_old, Wavefunction& wf_new, System & system, Kinetics &kinetics) {
    wf_old.update_old(wf_new, kinetics, moved_particle);
    system.update_old(moved_particle);
}

double Importance::get_g_ratio(const Wavefunction& wf_new, const Wavefunction & wf_old) {
    int j, i;
    double g_ratio;

    i = moved_particle;
    g_ratio = 0;
    for (j = 0; j < dim; j++) {
        g_ratio += 0.5 * (wf_old.qforce[i][j] + wf_new.qforce[i][j])*
                (D * timestep * 0.5 * (wf_old.qforce[i][j] - wf_new.qforce[i][j])
                - wf_new.r[i][j] + wf_old.r[i][j]);
    }

    return exp(g_ratio);
}

void Importance::set_trial_pos(Wavefunction& wf, System &system, Jastrow_factor & jastrow, Kinetics &kinetics) {
    int i, j;
    double qforce_test;

    qforce_test = 101;
    while (qforce_test > 100) {

        //Initial trial position
        for (i = 0; i < n_p; i++) {
            for (j = 0; j < dim; j++) {
                wf.r[i][j] = gaussian_deviate(&random_seed) * sqrt_timestep;
            }
        }

        wf.make_rel_matrix();

        system.initialize(wf, jastrow, kinetics);

        kinetics.get_QF(wf, system, jastrow);

        qforce_test = test_qforce(wf);

    }
}

double Importance::test_qforce(Wavefunction & wf) {
    int i, j;
    double test;

    test = 0;
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            if (abs(wf.qforce[i][j]) > test) {
                test = abs(wf.qforce[i][j]);
            }
        }
    }

    return test;
}

void Importance::set_new_pos(const Wavefunction &wf_old, Wavefunction & wf_new) {
    int i, j, k;

    i = moved_particle;

    //reseting old positions
    for (k = 0; k < n_p; k++) {
        if (k != i) {
            for (j = 0; j < dim; j++) {
                wf_new.r[k][j] = wf_old.r[k][j];
            }
        } else {
            for (j = 0; j < dim; j++) {
                wf_new.r[k][j] = wf_old.r[i][j] +
                        gaussian_deviate(&random_seed) * sqrt_timestep + wf_old.qforce[i][j] * timestep*D;
            }
        }
    }


    for (k = 0; k < n_p; k++) {
        for (j = k + 1; j < n_p; j++) {
            if ((j != i) && (k != i)) {
                wf_new.r_rel[k][j] = wf_new.r_rel[j][k] = wf_old.r_rel[k][j];
            }
        }
    }

    //updating the relative distance between the moved particle 
    //and the others
    for (j = 0; j < n_p; j++) {
        if (j != i) {
            wf_new.r_rel[i][j] = wf_new.r_rel[j][i] = wf_new.abs_relative(i, j);
        }
    }
}

void Importance::get_full_ratio(Wavefunction &wf_new, Wavefunction &wf_old, Jastrow_factor & jastrow) {
    j_ratio = jastrow.get_j_ratio(wf_new, wf_old, moved_particle);
    g_ratio = get_g_ratio(wf_new, wf_old);
    wf_ratio = s_ratio*j_ratio;

    full_ratio = g_ratio * wf_ratio*wf_ratio;
}

void Importance::do_vmc(Potential &pot, Kinetics &kinetics, System& system, Jastrow_factor& jastrow,
        Wavefunction& wf_old, Wavefunction & wf_new) {
    int cycles, i, j;
    
    //fjern
//    ofstream dist;
//    dist.open("dist_out.dat");

    jastrow.initialize_params(system);
    set_trial_pos(wf_old, system, jastrow, kinetics);

    for (cycles = 1; cycles <= n_c; cycles++) {
        for (moved_particle = 0; moved_particle < n_p; moved_particle++) {

            set_new_pos(wf_old, wf_new);
            s_ratio = system.get_spatial_ratio(wf_new, moved_particle);
         
            system.calc_for_newpos(wf_old, wf_new, kinetics, jastrow, moved_particle);
            kinetics.get_QF(wf_new, system, jastrow);

            get_full_ratio(wf_new, wf_old, jastrow);

            system.accepted_last = metropolis_test();
            if (system.accepted_last == true) {
                update_walker(wf_old, wf_new, system, kinetics);
//                if ((cycles > n_c/2) && (cycles%100 == 0)){
//                for (i=0; i < n_p; i++){
//                    for (j = 0; j < dim; j++){
//                        if (j == dim-1){
//                            dist << wf_old.r[i][j];
//                        } else{
//                                dist << wf_old.r[i][j] << " ";
//                        }
//                        }
//                    dist << endl;
//                }
//                }
            }

        }
        
        //if map_wf write_wf_map(cycles);

        delta_e = kinetics.get_KE(wf_old, system, jastrow);
        delta_e += pot.get_pot_E(wf_old);
        delta_e += pot.get_Coulomb(wf_old);

        local_e += delta_e;
        local_e2 += delta_e*delta_e;
    }
//    dist.close();

    local_e /= n_c;
    local_e2 /= n_c;

}

void Brute_Force::set_trial_pos(Wavefunction& wf, System& system, Jastrow_factor & jastrow, Kinetics &kinetics) {
    int i, j;
    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wf.r[i][j] = steplength * (ran1(&random_seed) - 0.5);
        }
    }
    wf.make_rel_matrix();

    if (kinetics.cf) {
        system.initialize(wf, jastrow, kinetics);
    }
    system.get_wf_val(wf, jastrow);
}

double fermions::get_det() {
    int i;
    double det, g1, g2;
    int *dummy = new int[n2];

    ludcmp(s_up, n2, dummy, &g1);
    ludcmp(s_down, n2, dummy, &g2);

    det = 1;
    for (i = 0; i < n2; i++) {
        det *= s_up[i][i] * s_down[i][i];
    }

    delete[] dummy;

    return g1 * g2*det;
}

double Pade_Jastrow::get_val(Wavefunction & wf) {
    int i, j;
    double arg;

    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            arg += a[i][j] * wf.r_rel[i][j] / (1.0 + beta * wf.r_rel[i][j]);
        }
    }

    return exp(arg);
}

void fermions::get_wf_val(Wavefunction& wf, Jastrow_factor & jastrow) {
    double value;

    initialize_slaters(wf);
    value = get_det();
    value *= jastrow.get_val(wf);

    wf.wf_val = value;
}

void Brute_Force::set_new_pos(const Wavefunction& wf_old, Wavefunction & wf_new) {
    int i, j;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wf_new.r[i][j] = wf_old.r[i][j] + steplength * (ran1(&random_seed) - 0.5);
        }
    }
    wf_new.make_rel_matrix();
}

void Brute_Force::get_full_ratio(Wavefunction& wf_new, Wavefunction& wf_old, Jastrow_factor & jastrow) {

    full_ratio = wf_new.wf_val / wf_old.wf_val;
    //full_ratio *= jastrow.get_j_ratio_full(wf_new, wf_old);
    full_ratio = full_ratio*full_ratio;
}

//double Pade_Jastrow::get_j_ratio_full(Wavefunction& wf_new, Wavefunction& wf_old) {
//    int i, j;
//    double j_ratio;
//
//    for (i = 0; i < n_p - 1; i++) {
//        for (j = i + 1; j < n_p; j++) {
//            j_ratio += a[i][j] * (wf_new.r_rel[i][j] / (1.0 + beta * wf_new.r_rel[i][j]) -
//                    wf_old.r_rel[i][j] / (1.0 + beta * wf_old.r_rel[i][j]));
//        }
//    }
//
//    return exp(j_ratio);
//}

void Brute_Force::update_walker(Wavefunction& wf_old, Wavefunction& wf_new, System & system, Kinetics &kinetics) {
    int i, j, k;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < dim; j++) {
            wf_old.r[i][j] = wf_new.r[i][j];
        }
        for (k = 0; k < n_p; k++) {

            wf_old.r_rel[k][i] = wf_old.r_rel[i][k] = wf_new.r_rel[k][i];
        }
    }

    wf_old.wf_val = wf_new.wf_val;
}

void Brute_Force::do_vmc(Potential &pot, Kinetics &kinetics, System& system, Jastrow_factor& jastrow, Wavefunction& wf_old, Wavefunction & wf_new) {
    int cycles;


    jastrow.initialize_params(system);

    set_trial_pos(wf_old, system, jastrow, kinetics);

    for (cycles = 1; cycles <= n_c; cycles++) {

        set_new_pos(wf_old, wf_new);
        system.get_wf_val(wf_new, jastrow);

        get_full_ratio(wf_new, wf_old, jastrow);

        system.accepted_last = metropolis_test();
        if (system.accepted_last == true) {
            update_walker(wf_old, wf_new, system, kinetics);

            //we need to initialize gradients and the inverse for the CF exp.
            if (kinetics.cf) {
                system.initialize(wf_old, jastrow, kinetics);
            }
        }


        delta_e = kinetics.get_KE(wf_old, system, jastrow);
        delta_e += pot.get_pot_E(wf_old);
        delta_e += pot.get_Coulomb(wf_old);

        local_e += delta_e;
        local_e2 += delta_e*delta_e;
    }

    local_e /= n_c;
    local_e2 /= n_c;
}

void initialize_qdots(int n_p, double w, int &dim, double &alpha, double &beta) {
    if (n_p == 2) {
        alpha = 0.987;
        beta = 0.398;
    } else if (n_p == 6) {
        if (w == 1) {
            alpha = 0.92;
            beta = 0.565;
        } else if (w == 0.28) {
            alpha = 0.88;
            beta = 0.33;
        } else {
            alpha = 0.64;
            beta = 0.09;
        }
    } else if (n_p == 12) {
        if (w == 1) {
            alpha = 0.87;
            beta = 0.68;
        } else {
            alpha = 0;
            beta = 0;
        }
    } else if (n_p == 20) {
        if (w == 1) {
            alpha = 0.84;
            beta = 0.76;
        } else {
            alpha = 0.87;
            beta = 0.68;
            cout << "warning" << endl;
        }

    } else if (n_p == 30) {
        if (w == 1) {
            alpha = 0.78;
            beta = 0.85;
        } else {

            alpha = 0.87;
            beta = 0.68;
            cout << "warning" << endl;
        }
    }
}

// random numbers with gaussian distribution

double gaussian_deviate(long * idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (idum < 0) iset = 0;
    if (iset == 0) {
        do {
            v1 = 2. * ran2(idum) - 1.0;
            v2 = 2. * ran2(idum) - 1.0;
            rsq = v1 * v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2. * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;

        return gset;
    }
} // end function for gaussian deviates

double get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec / 1000000;
    return d;
}

atoms::atoms(const VMC &vmc, const Nucleus &pot, double ALPHA) {
    n_p = vmc.n_p;
    n2 = vmc.n2;
    dim = vmc.dim;

    alpha = ALPHA;
    Z = pot.Z;

    a_sym = 0.25;
    a_asym = 0.5;

    inv_new = (double **) matrix(n2, n_p, sizeof (double));
    inv_old = (double **) matrix(n2, n_p, sizeof (double));

    s_down = (double **) matrix(n2, n2, sizeof (double));
    s_up = (double **) matrix(n2, n2, sizeof (double));
}

double atoms::phi(const Wavefunction& wf, int particle, int q_num) {
    double r, L, r_L, n;

    n = q_num + 1.0;
    r = sqrt(wf.get_r_i2(particle));
    r_L = 2 * alpha * r / n;

    //    x = wf.r[particle][0];
    //    y = wf.r[particle][1];


    if (q_num == 0) {
        L = 1;
    } else if (q_num == 1) {
        L = 2 - r_L;
    }//    } else if (q_num == 2) {
        //        L = 2 * y;
        //    }
        //    } else if (q_num == 3) {
        //        H = 4 * x * x - 2;
        //    } else if (q_num == 4) {
        //        H = 4 * y * y - 2;
        //    } else if (q_num == 5) {
        //        H = 4 * x*y;
        //    } else if (q_num == 6) {
        //        H = 4 * x * (2 * y * y - 1);
        //    } else if (q_num == 7) {
        //        H = 4 * y * (2 * x * x - 1);
        //    } else if (q_num == 8) {
        //        H = 8 * y * y * y - 12 * y;
        //    } else if (q_num == 9) {
        //        H = 8 * x * x * x - 12 * x;
        //    } else if (q_num == 10) {
        //        H = 8 * y * x * (2 * x * x - 3);
        //    } else if (q_num == 11) {
        //        H = 8 * y * x * (2 * y * y - 3);
        //    } else if (q_num == 12) {
        //        H = 16 * y * y * (y * y - 3) + 12;
        //    } else if (q_num == 13) {
        //        H = 16 * x * x * (x * x - 3) + 12;
        //    } else if (q_num == 14) {
        //        H = 4 * (2 * x * x - 1)*(2 * y * y - 1);
        //    } 
    else {
        cout << "Mismatching quantum number: " << q_num << endl;
    }

    return L * exp(-alpha * r / n);
}

double atoms::del_phi(const Wavefunction& wf, int particle, int q_num, int d) {

}

double atoms::lapl_phi(const Wavefunction& wf, int particle, int q_num) const {

}

Nucleus::Nucleus(const VMC& vmc, double z) {
    n_p = vmc.n_p;
    dim = vmc.dim;
    Z = z;
}

double Nucleus::get_pot_E(const Wavefunction& wf) const {
    int i;
    double r, e_potential;

    e_potential = 0;

    // contribution from oscillator part 
    for (i = 0; i < n_p; i++) {
        r = sqrt(wf.get_r_i2(i));
        e_potential -= 1 / r;
    }

    return Z*e_potential;
}

void initialize_atoms(int n_p, double Z, int& dim, double& alpha, double& beta) {
    dim = 3;
    alpha = 1;
    beta = 0.3;
}

//double LaguerreGeneral( int n, double alpha, double x)
//{
//  double *glaguerre; 
//  if ( alpha <= -1.0 ) {
//    cout << "LAGUERRE_GENERAL - Fatal error!" << endl;
//    cout << "The input value of ALPHA is=  "  << alpha  << endl;
//    cout << "but ALPHA must be greater than -1." << endl;
//  }
//  glaguerre = new double[n+1];
//  if ( n >= 0 ) {
//    for (int  i = 0; i <= n; i++) glaguerre[i] = 0.0;
//    glaguerre[0] = 1.0;
//    if ( n > 0 ){
//      glaguerre[1] = 1.0+alpha-x;
//      // recursion relation for generalized Laguerre polynomials
//      for (int  i = 2; i <= n; i++){
//          glaguerre[i] = ((2.0*i-1.0+alpha-x)*glaguerre[i-1]+
//	              (1.0-i-alpha)*glaguerre[i-2])/((float) i);
//      }
//    }
//  }
//  double GLaguerre = glaguerre[n];
//  delete [] glaguerre;
//  return GLaguerre;
//}   // end function glaguerre