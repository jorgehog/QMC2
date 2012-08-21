/* 
 * File:   vmc_oo.h
 * Author: jorgehog
 *
 * Created on 4. juli 2011, 14:08
 */

#ifndef VMC_OO_H
#define	VMC_OO_H

double gaussian_deviate(long* idum);
void initialize_qdots(int n_p, double w, int &dim,
        double &alpha, double &beta);
void initialize_atoms(int n_p, double Z, int &dim,
        double &alpha, double &beta);
double get_time();

class VMC;
class Kinetics;
class Potential;
class Harmonic_osc;
class Nucleus;
class Jastrow_factor;
class System;
class Wavefunction;

class VMC {
public:
    int n_c;

    //Wavefunction wf_old;
    //Wavefunction wf_new;
    
    int n_p;
    int n2;
    int dim;

    double delta_e;
    double local_e;
    double local_e2;

    double s_ratio;
    double j_ratio;
    double wf_ratio;
    double full_ratio;

    long random_seed;
    double random_number;
    
    bool map_wf;


    bool metropolis_test();
    void output();
    virtual void do_vmc(Potential &pot, Kinetics &kinetics, System &system, Jastrow_factor &jastrow,
            Wavefunction &wf_old, Wavefunction &wf_new) = 0;
    void free_memory(Wavefunction& wf_old, Wavefunction &wf_new, Jastrow_factor &jastrow, System &system, Kinetics &kinetics);
    virtual void set_trial_pos(Wavefunction &wf, System &system, Jastrow_factor &jastrow, Kinetics &kinetics) = 0;
    virtual void set_new_pos(const Wavefunction &wf_old, Wavefunction &wf_new) = 0;
    virtual void get_full_ratio(Wavefunction &wf_new, Wavefunction &wf_old, Jastrow_factor &jastrow) = 0;
    virtual void update_walker(Wavefunction &wf_old, Wavefunction &wf_new, System &system, Kinetics &kinetics) = 0;

};

class Importance : public VMC {
public:
    int moved_particle;
    double timestep;
    double sqrt_timestep;
    double D;
    double g_ratio;



    Importance(int N_C, int N_P, int DIM, double TIMESTEP);


    double test_qforce(Wavefunction &wf);
    double get_g_ratio(const Wavefunction &wf_new, const Wavefunction &wf_old);


    virtual void set_trial_pos(Wavefunction &wf, System &system, Jastrow_factor &jastrow, Kinetics &kinetics);

    virtual void do_vmc(Potential &pot, Kinetics &kinetics, System &system, Jastrow_factor &jastrow,
            Wavefunction &wf_old, Wavefunction &wf_new);
    virtual void set_new_pos(const Wavefunction &wf_old, Wavefunction &wf_new);

    virtual void get_full_ratio(Wavefunction &wf_new, Wavefunction &wf_old, Jastrow_factor &jastrow);
    virtual void update_walker(Wavefunction &wf_old, Wavefunction &wf_new, System &system, Kinetics &kinetics);



};

class Brute_Force : public VMC {
public:
    double steplength;

    Brute_Force(int N_C, int N_P, int DIM, double STEPLENGTH);

    virtual void set_trial_pos(Wavefunction &wf, System &system, Jastrow_factor &jastrow, Kinetics &kinetics);

    virtual void do_vmc(Potential &pot, Kinetics &kinetics, System &system, Jastrow_factor &jastrow,
            Wavefunction &wf_old, Wavefunction &wf_new);
    virtual void set_new_pos(const Wavefunction &wf_old, Wavefunction &wf_new);

    virtual void get_full_ratio(Wavefunction &wf_new, Wavefunction &wf_old, Jastrow_factor &jastrow);
    virtual void update_walker(Wavefunction &wf_old, Wavefunction &wf_new, System &system, Kinetics &kinetics);

};

class Jastrow_factor {
public:
    int n_p;
    int n2;
    int dim;

    double j_ratio;



    virtual double get_j_ratio(Wavefunction &wf_new, Wavefunction &wf_old, int i) = 0;
    //virtual double get_j_ratio_full(Wavefunction &wf_new, Wavefunction &wf_old)= 0;

    void get_grad(Wavefunction &wf);
    virtual double get_lapl_sum(const Wavefunction &wf) const = 0;
    virtual double get_deriv(Wavefunction &wf, int i, int k) = 0;
    virtual double get_val(Wavefunction &wf) = 0;
    virtual void initialize_params(const System &system) const = 0;

    virtual void clear() = 0;


};

class Pade_Jastrow : public Jastrow_factor {
public:
    double beta;
    double **a;

    Pade_Jastrow(const VMC &vmc, double BETA);


    virtual double get_j_ratio(Wavefunction &wf_new, Wavefunction &wf_old, int i);
    //virtual double get_j_ratio_full(Wavefunction &wf_new, Wavefunction &wf_old);

    virtual void initialize_params(const System &system) const;

    virtual double get_val(Wavefunction &wf);
    virtual double get_lapl_sum(const Wavefunction &wf) const;
    virtual double get_deriv(Wavefunction &wf, int i, int k);

    virtual void clear();
};

class System {
public:
    int n_p;
    int dim;
    double a_sym;
    double a_asym;
    double s_ratio;

    double** s_up;
    double** s_down;

    bool accepted_last;


    virtual double phi(const Wavefunction &wf, int particle, int q_num) = 0;
    virtual double del_phi(const Wavefunction &wf, int particle, int q_num, int d) = 0;
    virtual double lapl_phi(const Wavefunction &wf, int particle, int q_num) const = 0;
    virtual double get_spatial_ratio(const Wavefunction &wf_new, int i) = 0;

    virtual void initialize(Wavefunction &wf, Jastrow_factor &jastrow, Kinetics &kinetics) = 0;
    virtual void calc_for_newpos(const Wavefunction &wf_old, Wavefunction &wf_new, Kinetics &kinetics,
            Jastrow_factor &jastrow, int i) = 0;
    virtual void get_new_grad(const Wavefunction& wf_old, Wavefunction &wf_new, int particle) = 0;
    virtual void get_init_grad(Wavefunction &wf, Jastrow_factor &jastrow) = 0;
    virtual double get_lapl_sum(const Wavefunction &wf) const = 0;
    virtual void get_wf_val(Wavefunction &wf, Jastrow_factor &jastrow) = 0;
    virtual void update_old(int i) = 0;

    virtual void clear() = 0;


};

class fermions : public System {
public:
    int n2;

    double** inv_new;
    double** inv_old;

    double get_det();


    virtual double phi(const Wavefunction &wf, int particle, int q_num) = 0;
    virtual double del_phi(const Wavefunction &wf, int particle, int q_num, int d) = 0;
    virtual double lapl_phi(const Wavefunction &wf, int particle, int q_num) const = 0;
    virtual double get_spatial_ratio(const Wavefunction &wf_new, int i);

    void initialize_slaters(const Wavefunction &wf);
    void invert_slaters();
    void make_merged_inv(const Wavefunction &wf);
    virtual void initialize(Wavefunction &wf, Jastrow_factor &jastrow, Kinetics &kinetics);
    virtual void calc_for_newpos(const Wavefunction &wf_old, Wavefunction &wf_new, Kinetics &kinetics, Jastrow_factor &jastrow, int i);

    virtual void get_new_grad(const Wavefunction& wf_old, Wavefunction &wf_new, int particle);
    virtual void get_init_grad(Wavefunction &wf, Jastrow_factor &jastrow);
    virtual double get_lapl_sum(const Wavefunction &wf) const;
    virtual void get_wf_val(Wavefunction &wf, Jastrow_factor &jastrow);

    void update_old(int i);

    void clear();

};

class qdots : public fermions {
public:
    double alpha;
    double w;

    qdots(const VMC &vmc, const Harmonic_osc &pot, double ALPHA);

    virtual double phi(const Wavefunction &wf, int particle, int q_num);
    virtual double del_phi(const Wavefunction &wf, int particle, int q_num, int d);
    virtual double lapl_phi(const Wavefunction &wf, int particle, int q_num) const;



};

class atoms : public fermions {
public:
    double alpha;
    double Z;

    atoms(const VMC &vmc, const Nucleus &pot, double ALPHA);

    virtual double phi(const Wavefunction &wf, int particle, int q_num);
    virtual double del_phi(const Wavefunction &wf, int particle, int q_num, int d);
    virtual double lapl_phi(const Wavefunction &wf, int particle, int q_num) const;

};

class Wavefunction {
public:
    int n_p;
    int n2;
    int dim;

    double wf_val;

    double** r;
    double** r_rel;

    double** qforce;

    double** spatial_grad;
    double** jast_grad;


    //Wavefunction(const Wavefunction &wf);
    Wavefunction();
    Wavefunction(const VMC &vmc);


    double get_r_i2(int i) const;
    double abs_relative(int i, int j);

    void make_rel_matrix();
    void update_old(Wavefunction &wf_new, Kinetics &kinetics, int i);

    void clear();
};

class Kinetics {
public:
    int n_p;
    int dim;

    bool cf;

    virtual double get_KE(Wavefunction &wf, System &system, Jastrow_factor &jastrow) = 0;
    virtual void get_QF(Wavefunction &wf, System &system, Jastrow_factor &jastrow) = 0;
    virtual void clear() = 0;
};

class Kinetics_num : public Kinetics {
public:
    double h, h2;

    Wavefunction wfplus;
    Wavefunction wfminus;

    Kinetics_num(VMC &vmc);
    Kinetics_num(VMC &vmc, double H);

    virtual double get_KE(Wavefunction &wf, System &system, Jastrow_factor &jastrow);
    virtual void get_QF(Wavefunction &wf, System &system, Jastrow_factor &jastrow);

    virtual void clear();
};

class Kinetics_cf : public Kinetics {
public:
    Kinetics_cf(VMC &vmc);

    virtual double get_KE(Wavefunction &wf, System &system, Jastrow_factor &jastrow);
    virtual void get_QF(Wavefunction &wf, System &system, Jastrow_factor &jastrow);

    virtual void clear();

};

//class Energy {
//public:
//    int n_p;
//    int dim;
//
//    double e_potential;
//    double e_kinetic;
//
//
//    double h, h2;
//
//    Energy();
//    Energy(const VMC &vmc);
//
//    virtual double get_pot_E(const Wavefunction &wf) const = 0;
//    double get_LE(const Wavefunction &wf, const System &system, const Jastrow_factor &jastrow);
//    double get_LE_num(const Wavefunction& wf, System& system, Jastrow_factor &jastrow);
//
//};

class Potential {
public:
int n_p;
int dim;

virtual double get_pot_E(const Wavefunction &wf) const = 0;
double get_Coulomb(const Wavefunction &wf) const;
};

class Harmonic_osc : public Potential {
public:
    double w;

    Harmonic_osc(const VMC &vmc, double W);

    virtual double get_pot_E(const Wavefunction& wf) const;

};

class Nucleus : public Potential {
public:
    double Z;

    Nucleus(const VMC &vmc, double z);

    virtual double get_pot_E(const Wavefunction& wf) const;

};

#endif	/* VMC_OO_H */



