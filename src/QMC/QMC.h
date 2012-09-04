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

    int cycle;

    int accepted;
    int thermalization;

    double local_E;

    Walker* original_walker;
    Walker* trial_walker;

    Jastrow *jastrow;
    Sampling *sampling;
    System *system;
    Kinetics *kinetics;

    std::vector<OutputHandler*> output_handler;

    QMC(int n_p, int dim, int n_c,
            Sampling *sampling,
            System *system,
            Kinetics *kinetics = new NoKinetics(),
            Jastrow *jastrow = new No_Jastrow());


    virtual void initialize() = 0;

    void dump_output();
    void finalize_output();

    void diffuse_walker();

    void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const;
    double get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const;

    void calculate_energy_necessities(Walker* walker) const;

    bool metropolis_test(double A);
    void update_walker(Walker*walker_pre, const Walker* walker_post, int particle) const;
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

    void copy_walker(const Walker* parent, Walker* child) const;

public:

    void add_output(OutputHandler* output_handler);

    virtual void run_method() = 0;
    virtual void output() const = 0;

    void get_gradients(Walker* walker, int particle) const;
    void get_gradients(Walker* walker) const;
    void get_wf_value(Walker* walker) const;
    void get_laplsum(Walker* walker) const;

    void update_pos(const Walker* walker_pre, Walker* walker_post, int particle) const;
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
        return accepted / double(n_p * (n_c + thermalization));
    }
    
    friend class Distribution;

};

class VMC : public QMC {
protected:

    double vmc_E, E2;

    virtual void initialize();
    void calculate_energy(Walker* walker);
    void scale_values();

public:

    VMC(int n_p, int dim, int n_c,
            Sampling *sampling,
            System *system,
            Kinetics *kinetics = new NoKinetics,
            Jastrow *jastrow = new No_Jastrow()
            );

    double get_var() const;
    double get_energy() const;
    double get_e2() const;
    void set_e(double e);
    void set_e2(double e2);

    virtual void run_method();
    virtual void output() const;

    friend class Minimizer;
    friend class ASGD;

    friend class BlockingData;

};

class DMC : public QMC {
protected:

    int K; //Factor of empty space for walkers over initial walkers
    int n_w_orig;
    int n_w;
    int n_w_last;

    int block_size;
    int samples;

    double dmc_E;
    double E_T;
    double E;

    bool dist_from_file;

    Walker **Angry_mob;

    void initialize();

    void initialize_walker(int k);
    void Evolve_walker(double GB);

    void bury_the_dead();

    void update_energies();

public:

    DMC(int n_p, int dim, int n_w, int n_c, int block_size, double E_T,
            Sampling *sampling,
            System *system,
            Kinetics *kinetics = new NoKinetics(),
            Jastrow *jastrow = new No_Jastrow(),
            bool dist_from_file = false);

    virtual void run_method();

    virtual void output() const;
    
    friend class stdoutDMC;

};



#endif	/* QMC_H */
