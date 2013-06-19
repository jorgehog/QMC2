/* 
 * File:   Fermions.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:44 PM
 */

#ifndef FERMIONS_H
#define	FERMIONS_H

#include "../System.h"

#include <armadillo>
#include "../../Walker/Walker.h"
struct GeneralParams;


/*! \brief The Fermion system class. */
class Fermions : public System {
protected:

    //! The diagonal of the new slater matrix times the old slater inverse.
    /*! 
     * Needed for updating the inverse.
     * Stored because only half of the vector is changed when moving one particle.
     * \see System::set_spin_state()
     */
    arma::rowvec I;
    bool node_crossed; //!< True if the spatial ratio is negative.

    //! Method for calculating the Slater matrix inverse
    /*!
     * The merged inverse is made by concatenating the two slater matrix inverses.  
     * This way we can sum freely over particles without having to if-test on the spin.
     */
    void make_merged_inv(Walker* walker);

    //!Method for updating the inverse given that we moved one particle.
    void update_inverse(const Walker* walker_old, Walker* walker_new, int particle);


public:

    Fermions(GeneralParams &, Orbitals*);

    void get_spatial_grad(Walker* walker, int particle) const;
    void get_spatial_grad_full(Walker* walker) const;
    double get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle);
    double get_spatial_lapl_sum(Walker* walker) const;

    //! Fixed node approximation.
    bool allow_transition() {
        return !node_crossed;
    }

    /*!
     * Copies the inverse.
     */
    void copy_walker(const Walker* parent, Walker* child) const {
        child->inv = parent->inv;
    }

    /*!
     * Updates the inverse.
     */
    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
        using namespace arma;
        walker_pre->inv(span(0, n2 - 1), span(start, end)) = walker_post->inv(span(0, n2 - 1), span(start, end));
    };

    /*!
     * Resets the inverse.
     */
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
        using namespace arma;
        walker_post->inv(span(0, n2 - 1), span(start, end)) = walker_pre->inv(span(0, n2 - 1), span(start, end));
    }

    /*!
     * The determinant of each spin value multiplied.
     */
    double get_spatial_wf(const Walker* walker) {
        using namespace arma;
        return det(walker->phi(span(0, n2 - 1), span())) * det(walker->phi(span(n2, n_p - 1), span()));
    }

    /*!
     * Calculates the inverse.
     */
    void initialize(Walker* walker) {
        make_merged_inv(walker);
    }

    /*!
     * When a particle is moved, the inverse is updated.
     */
    void calc_for_newpos(const Walker* walker_old, Walker* walker_new, int i) {
        update_inverse(walker_old, walker_new, i);
    }

};

#endif	/* FERMIONS_H */

