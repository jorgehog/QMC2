#pragma once


#include "../Sampling.h"

/*! \brief Implementation of the Brute Force QMC.
 * Uses the Simle diffusion class.
 * All methods are empty except for the energy necessities part which requires the
 * gradients to be calculated (not using the Quantum Force).
 */
class Brute_Force : public Sampling {
public:

    Brute_Force(GeneralParams &);

    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
        (void) walker_pre;
        (void) walker_post;
        (void) particle;
    }

    void get_necessities(Walker* walker){
        (void) walker;
    }
    
    void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
        (void) walker_pre;
        (void) walker_post;
        (void) particle;
    }

    void calculate_energy_necessities(Walker* walker) const {
        qmc->get_gradients(walker);
    }

    void copy_walker(const Walker* parent, Walker* child) const {
        (void) parent;
        (void) child;
    }
    
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
        (void) walker_pre;
        (void) walker_post;
        (void) particle;
    }

};

