/* 
 * File:   stdoutASGD.h
 * Author: jorgmeister
 *
 * Created on October 29, 2012, 3:05 PM
 */

#ifndef STDOUTASGD_H
#define	STDOUTASGD_H

class stdoutASGD : public OutputHandler {
public:
    stdoutASGD(std::string filename = "ASGD_out",
            std::string path = "./");

    virtual void dump();
    
    virtual void post_pointer_init(){
        grad = arma::zeros(asgd->Nparams);
    }
    
private:
    arma::vec grad;
    double sumE;
};

#endif	/* STDOUTASGD_H */

