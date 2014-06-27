#pragma once


#include <armadillo>
#include <string>


namespace QMC2
{


/*! \brief Class representing a Random Walker.
 * Holds position data, alive state, etc.
 * Designed to lighten function arguments, and ease implementation
 * of QMC methods involving multiple walkers.
 * Alot of values are stored to avoid calculating the same value twice.
 */
class Walker {
protected:
    int n_p;
    int n2;
    int dim;

    bool is_murdered; //!< If true, the walker will be deleted and removed (DMC only).

public:

    const int &_n_p() const
    {
        return n_p;
    }

    const int &_dim() const
    {
        return dim;
    }


    //! Constructor.
    /*!
     * @param alive If false, the walker is initialized dead.
     */
    Walker(int n_p, int dim, bool alive = true);

    double spatial_ratio; //!< The ratio of the spatial wave function (stored in the newest walker).
    double lapl_sum; //!< The sum of the Laplacians of all particles.
    double E; //!< The energy of the given configuration (stored to speed up DMC).

    arma::mat r; //!< The positions of all particles.
    arma::mat r_rel; //!< The relative positions of all particles.

    arma::mat qforce; //!< The Quantum Force for all particles.

    arma::mat spatial_grad; //!< The gradient of the Spatial Wave function for all particles.
    arma::mat jast_grad; //!< The gradient of the Jastrow Factor for all particles.
    arma::mat inv; //!< The inverse of the Slater matrix (given fermion system)

    arma::mat phi; //!< The single particle wave functions for all particles and quantum numbers.
    arma::field<arma::mat> dell_phi; //!< The derivatives of the single particle wave functions for all particles and quantum numbers.
    arma::cube dJ; //!< Cube used for storing sum terms for the Jastrow Factor's closed form expressions. 

    arma::rowvec r2; //!< The radius squared for all particles.
    arma::rowvec abs_r; //!< The radius for all particles;


    //! Method for calculating the radius squared for one particle.
    /*!
     * @param i The particle number.
     */
    void calc_r_i2(int i);

    //! Method for calculating the radius squared for all particles.

    void calc_r_i2() {
        for (int i = 0; i < n_p; i++) {
            this->calc_r_i2(i);
        }
    }

    //! Method for calculating the radius of a particle. Assumes the squared exist.
    /*!
     * @param i Particle number.
     */
    void calc_r_i(int i) {
        abs_r(i) = sqrt(get_r_i2(i));
    };

    //! Method for calculating the radius for all particles.

    void calc_r_i() {
        for (int i = 0; i < n_p; i++) {
            this->calc_r_i(i);
        }
    }

    //! Method for calculating the relative distance between two particles.
    /*!
     * @param i,j The particle numbers.
     */
    double calc_r_rel(int i, int j) const;

    /*!
     * Creates the relative position matrix.
     */
    void make_rel_matrix();

    /*!
     * Send a walker to a different node.
     * @param dest The receiving node's rank. 
     */
    void send_soul(int dest);

    /*!
     * Receives a walker from a different node.
     * @param root The rank of the node from which the walker was sent.
     */
    void recv_soul(int root);

    //! Method for fetching the squared radius of a particle. 

    /*!
     * Used in order to avoid calculating the same radius twice.
     * @param i Particle number.
     */
    double get_r_i2(int i) const {
        return r2(i);
    }

    //! Method for calculating the radius of a particle.

    /*!
     * @param i Particle number.
     */
    double get_r_i(int i) const {
        return abs_r(i);
    }

    /*!
     * Flags the walker for destruction.
     * \see DMC::bury_the_dead()
     */
    void kill() {
        is_murdered = true;
    }

    bool is_dead() {
        return is_murdered;
    }

    bool is_alive() {
        return !is_murdered;
    }

    /*!
     * Sets the destruction flag to false.
     */
    void ressurect() {
        is_murdered = false;
    }

    void set_E(double E) {
        this->E = E;
    }

    double get_E() const {
        return E;
    }

    //! Prints out all the walkers information.
    /*!
     * Extremely handy for debugging.
     * @param header A header for the printout in order to distinguish several printouts easily.
     */
    void print(std::string header = "----") const;


};

}
