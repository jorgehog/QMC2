/* 
 * File:   Bosons.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:44 PM
 */

#ifndef Bosons_H
#define	Bosons_H

/*! \brief The Boson system class. */
class Bosons : public System {
protected:
    
    int a; //!< The hard core radius for the infinite potential.
    bool overlap; //!< True if the relative distance is less than the hard core radius.

public:

    Bosons(GeneralParams &, Orbitals*);

    void get_spatial_grad(Walker* walker, int particle) const;
    void get_spatial_grad_full(Walker* walker) const;
    double get_spatial_ratio(const Walker* walker_post, const Walker* walker_pre, int particle);
    double get_spatial_lapl_sum(const Walker* walker) const;

    //! Infinite potential to simulate bosonic behaviour.
    bool allow_transition() {
        return !overlap;
    }

    /*!
     * Does nothing.
     */
    void copy_walker(const Walker* parent, Walker* child) const {

    }

    /*!
     * Does nothing.
     */
    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {

    };

    /*!
     * Does nothing.
     */
    void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
    
    }

    /*!
     * The single particle states of each particle multiplied.
     * Assumes the trial wave function initializes every particle in the same
     * single particle state.
     */
    double get_spatial_wf(const Walker* walker);

    /*!
     * Does nothing.
     */
    void initialize(Walker* walker) {
        
    }

    /*!
     * Does nothing.
     */
    void calc_for_newpos(const Walker* walker_old, Walker* walker_new, int i) {

    }

};

#endif	/* Bosons_H */

