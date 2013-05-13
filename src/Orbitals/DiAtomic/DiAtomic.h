/* 
 * File:   DiAtomic.h
 * Author: jorgmeister
 *
 * Created on May 10, 2013, 1:46 PM
 */

#ifndef DIATOMIC_H
#define	DIATOMIC_H

class DiAtomic : public Orbitals {
public:

    DiAtomic(GeneralParams & gP, VariationalParams & vP);

    /*!
     * Calculates the exponential terms exp(-r/n) for all needed n once pr. particle
     * per core to save CPU-time. 
     * \see Orbitals::set_qnum_indie_terms()
     */
    void set_qnum_indie_terms(Walker * walker, int i);

protected:

    double* R;

    hydrogenicOrbitals* nucleus1;
    hydrogenicOrbitals* nucleus2;

    Walker* walker_nucleus1;
    Walker* walker_nucleus2;

    double getAlpha(VariationalParams & vP, double k0);
    
    ///! Sums contrib from nucleus1 and 2 in their mass center coordinates.
    double get_variational_derivative(const Walker* walker, int n);

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

#endif	/* DIATOMIC_H */

