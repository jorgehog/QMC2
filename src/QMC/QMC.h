/* 
 * File:   QMC.h
 * Author: jorgehog
 *
 * Created on 30. mars 2012, 17:42
 */

#ifndef QMC_H
#define	QMC_H

class QMC {
protected:
    int n_c;

    int n_p;
    int n2;
    int dim;

    int accepted;

    Walker* original_walker;
    Walker* trial_walker;

    Jastrow *jastrow;
    Sampling *sampling;
    System *system;
    Kinetics *kinetics;

    QMC(int n_p, int dim, int n_c, Jastrow *jastrow, Sampling *sampling, System *system, Kinetics *kinetics);


    virtual void initialize() = 0;

    void diffuse_walker();

    void update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const;
    void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const;
    double get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const;

    void calculate_energy_necessities(Walker* walker) const;

    bool metropolis_test(double A);
    void update_walker(Walker*walker_pre, const Walker* walker_post, int particle) const;
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

    void copy_walker(const Walker* parent, Walker* child) const;

public:

    virtual void run_method() = 0;
    virtual void output() const = 0;

    void get_gradients(Walker* walker, int particle) const;
    void get_gradients(Walker* walker) const;
    void get_wf_value(Walker* walker) const;
    void get_laplsum(Walker* walker) const;

    double calculate_local_energy(Walker* walker) const;

    System* get_system_ptr() {
        return system;
    }

    Kinetics* get_kinetics_ptr() {
        return kinetics;
    }

    Sampling* get_sampling_ptr() {
        return sampling;
    }

    Jastrow* get_jastrow_ptr() {
        return jastrow;
    }

    Orbitals* get_orbitals_ptr() {
        return system->get_orbital_ptr();
    }

    double get_accepted_ratio() const {
        return accepted / double(n_p);
    }

    friend class Minimizer;
    friend class ASGD;

};

class VMC : public QMC {
protected:
    double vmc_E, E2;

    bool dist_to_file;

    virtual void initialize();
    void calculate_energy(Walker* walker);
    void scale_values();

public:

    VMC(int n_p, int dim, int n_c, Jastrow *jastrow, Sampling *sampling,
            System *system, Kinetics *kinetics, bool dist_to_file = false);

    double get_var() const;
    double get_energy() const;
    double get_e2() const;
    void set_e(double e);
    void set_e2(double e2);

    virtual void run_method();
    virtual void output() const;

    friend class Minimizer;
    friend class ASGD;

};

class DMC : public QMC {
protected:

    int K; //Factor of empty space for walkers over initial walkers
    int n_w_orig, n_w, n_w_last;
    int samples;
    double E_T, E, E_avg;
    bool dist_from_file;

    Walker **Angry_mob;

    virtual void initialize();
    void initialize_walker(int k);
    void increase_walker_space();
    void bury_the_dead();
    void Evolve_walker(double GB);
    void update_energies(int n);

    //problem: DMC not equal for IMPORTANCE vs. BF
public:
    DMC(int n_p, int dim, int n_w, int n_c, double E_T, Jastrow *jastrow, Sampling *sampling,
            System *system, Kinetics *kinetics, bool dist_from_file = false);


    virtual void run_method();

    virtual void output() const;
};



#endif	/* QMC_H */
