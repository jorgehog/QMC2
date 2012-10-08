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

    Jastrow *jastrow;
    Sampling *sampling;
    System *system;
    Kinetics *kinetics;

    std::vector<OutputHandler*> output_handler;

    virtual void initialize() = 0;

    virtual bool move_autherized(double A) = 0;

    void dump_output();
    void finalize_output();

    void diffuse_walker(Walker* original, Walker* trial);

    double get_acceptance_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const;

    void calculate_energy_necessities(Walker* walker) const;

    bool metropolis_test(double A);
    void update_walker(Walker*walker_pre, const Walker* walker_post, int particle) const;
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

    void copy_walker(const Walker* parent, Walker* child) const;

public:

    QMC(int n_p, int dim, int n_c,
            Sampling *sampling,
            System *system,
            Kinetics *kinetics = new NoKinetics(),
            Jastrow *jastrow = new No_Jastrow());
    QMC();

    void add_output(OutputHandler* output_handler);

    virtual void run_method() = 0;
    virtual void user_output() const = 0;

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

    double get_accepted_ratio(int total_cycles) const {
        return accepted / double(total_cycles);
    }


};

class VMC : public QMC {
protected:

    double vmc_E, E2;

    Walker* original_walker;
    Walker* trial_walker;

    virtual void initialize();

    virtual bool move_autherized(double A);

    void calculate_energy(Walker* walker);
    void scale_values();

public:

    VMC(GeneralParams &, VMCparams &, SystemObjects &);

    double get_var() const;
    double get_energy() const;
    double get_e2() const;
    void set_e(double e);
    void set_e2(double e2);

    virtual void run_method();
    virtual void user_output() const;

    friend class Minimizer;
    friend class ASGD;

    friend class Distribution;
    friend class BlockingData;

};

class DMC : public QMC {
protected:

    int K; //Factor of empty space for walkers over initial walkers

    int n_w;
    int n_w_last;

    int deaths;

    int block_size;
    int samples;

    double dmc_E;
    double E_T;
    double E;

    bool dist_from_file;
    std::string dist_in_path;

    Walker **original_walkers;
    Walker *trial_walker;

    void initialize();

    virtual bool move_autherized(double A);

    void iterate_walker(int k, int n_b = 1);

    void Evolve_walker(int k, double GB);

    void bury_the_dead();

    void update_energies();

    void reset_parameters() {
        n_w_last = n_w;
        E = samples = deaths = 0;
    }

public:

    DMC(GeneralParams &, DMCparams &, SystemObjects &);

    virtual void run_method();

    virtual void user_output() const;

    friend class stdoutDMC;

};



#endif	/* QMC_H */
