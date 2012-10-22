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
            Jastrow *jastrow = new No_Jastrow());
    QMC();

    void add_output(OutputHandler* output_handler);

    virtual void run_method() = 0;
    virtual void user_output() const = 0;

    double get_KE(const Walker* walker) const;
    void get_QF(Walker* walker) const;

    void get_gradients(const Walker* walker_pre, Walker* walker_post, int particle) const;
    void get_gradients(Walker* walker) const;
    void get_wf_value(Walker* walker) const;
    void get_laplsum(Walker* walker) const;

    double calculate_local_energy(const Walker* walker) const {
        return get_KE(walker) + system->get_potential_energy(walker);
    }

    System* get_system_ptr() const {
        return system;
    }

    Sampling* get_sampling_ptr() const {
        return sampling;
    }

    Jastrow* get_jastrow_ptr() const {
        return jastrow;
    }

    Orbitals* get_orbitals_ptr() const {
        return system->get_orbital_ptr();
    }

    double get_accepted_ratio(int total_cycles) const {
        return accepted / double(total_cycles);
    }

};

#endif	/* QMC_H */
