#pragma once


#include "../OutputHandler.h"

class DMC;

/*! \brief Class for handling the output of DMC.
 * Outputs values such as the trial energy, dmc energy, number of walkers, etc.
 */
class stdoutDMC : public OutputHandler {
protected:

    int n; //!< Number of times the dump() method has been called.
    double sumE; //!< Sum of the DMC energy used to calculate the trailing average.
    double sumN; //!< Sum of the number of walkers used to calculate the trailing average.

    DMC* dmc;
    
public:
    
    stdoutDMC(DMC* dmc, std::string path);

    void dump();

    void reset();

};

