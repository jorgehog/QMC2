/* 
 * File:   stdoutASGD.h
 * Author: jorgmeister
 *
 * Created on October 29, 2012, 3:05 PM
 */

#ifndef STDOUTASGD_H
#define	STDOUTASGD_H

/*!
 * \brief Class for handling the output of ASGD. 
 * Ouputs values such as the variational gradients, step length, variational parameters, etc.
 */
class stdoutASGD : public OutputHandler {
public:
    
    stdoutASGD(std::string path, std::string filename = "ASGD_out");

    void dump();
    
    /*!
     * Initializes the correct size of the variational gradient once the min
     * pointer has been cast to ASGD.
     */
    void post_pointer_init(){
        grad = arma::zeros(asgd->Nparams);
    }
    
private:
    
    arma::vec grad; //!< Vector used for calulation the trailing averages of variational parameters.
    double sumE; //!< Sum of the sampled energy used to calculate the trailing average.

};

#endif	/* STDOUTASGD_H */

