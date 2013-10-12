/* 
 * File:   DiTransform.h
 * Author: jorgmeister
 *
 * Created on May 10, 2013, 1:46 PM
 */

#ifndef DITRANSFORM_H
#define	DITRANSFORM_H

#include "../Orbitals.h"
struct GeneralParams;
struct VariationalParams;


class DiTransform : public Orbitals {
public:

    DiTransform(GeneralParams & gP, VariationalParams & vP, int system);

    /*!
     * Calculates the exponential terms exp(-r/n) for all needed n once pr. particle
     * per core to save CPU-time. 
     * \see Orbitals::set_qnum_indie_terms()
     */
    void set_qnum_indie_terms(Walker * walker, int i);

    void debug();

protected:

    double* R;

    Orbitals* nucleus1;
    Orbitals* nucleus2;

    Walker* walker_nucleus1;
    Walker* walker_nucleus2;

    double static minusPower(int n) {
        return -2 * (n % 2) + 1;
    }

    ///! Sums contrib from nucleus1 and 2 in their mass center coordinates.
    double get_dell_alpha_phi(Walker* walker, int p, int q_num);

    /*!
     * @return The variational parameter alpha for all objects.
     */
    double get_parameter(int n) {
        return nucleus1->get_parameter(n);
    }

    /*!
     * Calls methods in hydrogenicOrbitals.
     */
    void set_parameter(double parameter, int n) {
        nucleus1->set_parameter(parameter, n);
        nucleus2->set_parameter(parameter, n);
    }

    double phi(const Walker* walker, int particle, int q_num);

    //! Calculates the single particle wave function derivative for a given walker's particle and dimension.
    /*!
     * @param q_num The quantum number index.
     * @param d The dimension for which the derivative should be calculated (x,y,z).
     */
    double del_phi(const Walker* walker, int particle, int q_num, int d);

    //! Calculates the single particle wave function for a given walker's particle.
    /*!
     * @param q_num The quantum number index.
     */
    double lapl_phi(const Walker* walker, int particle, int q_num);

};

#endif	/* DITRANSFORM_H */

