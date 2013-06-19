/* 
 * File:   Importance.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:43 PM
 */

#ifndef IMPORTANCE_H
#define	IMPORTANCE_H

#include "../Sampling.h"

/*! \brief Implementation of Importance sampled QMC.
 * Using the Fokker-Planck diffusion class.
 * Introduces the Quantum Force.
 */
class Importance : public Sampling {
public:

    Importance(GeneralParams &);

    /*!
     * The parts of the gradients with the same spin as the moved particle are updated.
     */
    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const;

    /*!
     * The parts of the gradients with the same spin as the moved particle are reset.
     */
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const;

    /*!
     * The gradients and the Quantum force are calculated.
     */
    void get_necessities(Walker* walker) {
        qmc->get_gradients(walker);
        qmc->get_QF(walker);
    }

    /*!
     * The gradients are updated and the Quantum force is re-calculated.
     */
    void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
        qmc->get_gradients(walker_pre, walker_post, particle);
        qmc->get_QF(walker_post);
    }

    /*!
     * No energy necessities (they are already calculated).
     */
    void calculate_energy_necessities(Walker* walker) const {

    }

    /*!
     * The gradients and the Quantum force is copied.
     */
    void copy_walker(const Walker* parent, Walker* child) const;


};

#endif	/* IMPORTANCE_H */

