/* 
 * File:   SimpleVar.h
 * Author: jorgmeister
 *
 * Created on October 30, 2012, 6:13 PM
 */

#ifndef SIMPLEVAR_H
#define	SIMPLEVAR_H

class SimpleVar : public ErrorEstimator {
public:
    SimpleVar(int n_c, ParParams &);
    
    double estimate_error(){
        std::cout << mean(data) << std::endl;
        return arma::var(data);
    }
    
protected:
    
    double E;
    double E2;
    
    
};

#endif	/* SIMPLEVAR_H */

