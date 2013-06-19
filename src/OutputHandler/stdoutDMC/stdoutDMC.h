/* 
 * File:   stdoutDMC.h
 * Author: jorgmeister
 *
 * Created on September 4, 2012, 7:22 PM
 */

#ifndef STDOUTDMC_H
#define	STDOUTDMC_H

#include "../OutputHandler.h"

/*! \brief Class for handling the output of DMC.
 * Outputs values such as the trial energy, dmc energy, number of walkers, etc.
 */
class stdoutDMC : public OutputHandler {
protected:

    int n; //!< Number of times the dump() method has been called.
    double sumE; //!< Sum of the DMC energy used to calculate the trailing average.
    double sumN; //!< Sum of the number of walkers used to calculate the trailing average.

public:

    stdoutDMC(std::string path, std::string filename = "DMC_out");

    void dump();

};

#endif	/* STDOUTDMC_H */

