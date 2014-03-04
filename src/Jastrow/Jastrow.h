#pragma once


namespace QMC2
{


class Walker;

/*! \brief The class representing the Jastrow correlation functions
 * Holds all data concerning the Jastrow function and it's influence
 * on the QMC algorithm. 
 */
class Jastrow {
protected:
    int n_p;
    int n2;
    int dim;

    bool active; //!< Parameter false if No_Jastrow is loaded.

    //! Returns variational parameters
    /*!
     * @param n The index of the sought variational parameter
     * @return Variational parameter with index [n]
     */
    virtual double get_parameter(int n) = 0;

    //! Sets variational parameters
    /*!
     * @param n The index of the sought variational parameter
     * @param param The new value of parameter [n]
     */
    virtual void set_parameter(double param, int n) = 0;

    //! Calculates the derivative of the Jastrow factor with respect to a variational parameter
    /*!
     * @param n The index of the variational parameter for which the derivative is to be taken
     * @param walker The walker holds the positions etc. needed to evaluate the derivative 
     */
    virtual double get_variational_derivative(const Walker* walker, int n);

    //! Numerical Cartesian derivative. 
    /*!
     * For use in get_grad() when no closed form expression is implemented.
     * @param i Particle number.
     * @param d Dimension (x,y,z).
     */
    double get_derivative_num(Walker* walker, int i, int d) const;

    //! Numerical Cartesian Laplacian.
    /*!
     * For use in get_lapl_sum() when no closed form expression is implemented.
     */
    double get_laplaciansum_num(Walker* walker) const;

public:
    Jastrow(int n_p, int dim);
    Jastrow();

    /*!
     * Initializes the non-variational parameters needed by the Jastrow Factor.
     */
    virtual void initialize() = 0;

    /*!
     * Calculates the value of the Jastrow Factor at the walker's position.
     */
    virtual double get_val(const Walker* walker) const = 0;

    //! Calculates the ratio of the Jastrow factor needed by metropolis.
    /*!
     * @param walker_new Walker at current time step
     * @param walker_old Walker at previous time step
     * @param i The particle number.
     */
    virtual double get_j_ratio(const Walker* walker_new, const Walker* walker_old, int i) const = 0;

    /*!
     * Calculates the entire Cartesian gradient.
     */
    virtual void get_grad(Walker* walker) const = 0;

    //! Updates the gradient for a new particle move.
    /*!
     * @param walker_post Walker at current time step
     * @param walker_pre Walker at previous time step
     * @param i Particle number.
     */
    virtual void get_grad(const Walker* walker_pre, Walker* walker_post, int i) const = 0;

    //! Updates the summation factors of Jastrow factors closed form expressions.
    /*!
     * Used to optimize the calculations as few of these terms change as we move a particle. 
     * @param i Particle number.
     */
    virtual void get_dJ_matrix(Walker* walker, int i) const = 0;

    //! Calculates the summation factors of Jastrow factors closed form expressions.
    /*!
     * Used to optimize the calculations as few of these terms change as we move a particle. 
     */
    void get_dJ_matrix(Walker* walker) const;

    //! Method for calculating the Laplacian.
    /*!
     * Calculates the sum of all particles Laplacians.
     */
    virtual double get_lapl_sum(Walker* walker) const = 0;


    friend class Minimizer;
    friend class ASGD;
    friend class stdoutASGD;

};

}
