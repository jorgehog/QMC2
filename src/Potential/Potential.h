#pragma once


#include <string>
#include "../Sampler/Sampler.h"


namespace QMC2
{


class Walker;

/*! \brief Superclass for potentials. 
 * Potentials are stores in a vector in the system object.
 * \see System::potentials, System::get_potential_energy()
 */
class Potential {
protected:

    int n_p;
    int dim;

    std::string name;

public:


    Potential(int n_p, int dim);
    Potential();

    //! Method for calculating a walker's potential energy.
    /*!
     * Method overridden by subclasses.
     * @param walker The walker for which the potential energy should be calculated.
     * @return The potential energy.
     */
    virtual double get_pot_E(const Walker* walker) const = 0;

    Sampler pot_sampler;
    
    std::string get_name(){
        return name;
    }

};

}
