#pragma once


#include "../OutputHandler.h"

#include <armadillo>


namespace QMC2
{

class ASGD;

/*!
 * \brief Class for handling the output of ASGD. 
 * Outputs values such as the variational gradients, step length, variational parameters, etc.
 */
class stdoutASGD : public OutputHandler {
public:
    
    stdoutASGD(ASGD* asgd, std::string path);

    void dump();

    void reset();
    
    
private:
    
    ASGD* asgd;
    arma::vec grad; //!< Vector used for calulation the trailing averages of variational parameters.
    double sumE; //!< Sum of the sampled energy used to calculate the trailing average.

};

}
