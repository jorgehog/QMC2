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
    SimpleVar(int n_c);
    
    double estimate_error();
    
protected:
    
    double E;
    double E2;
    
    
};

#endif	/* SIMPLEVAR_H */

