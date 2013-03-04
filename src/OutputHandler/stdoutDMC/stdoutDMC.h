/* 
 * File:   stdoutDMC.h
 * Author: jorgmeister
 *
 * Created on September 4, 2012, 7:22 PM
 */

#ifndef STDOUTDMC_H
#define	STDOUTDMC_H

class stdoutDMC : public OutputHandler {
protected:

    int n;
    double sumE;
    double sumN;

public:

    stdoutDMC(std::string path, std::string filename = "DMC_out");

    virtual void dump();

};

#endif	/* STDOUTDMC_H */

