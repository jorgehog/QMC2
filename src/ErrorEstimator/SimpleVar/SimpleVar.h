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
    SimpleVar(ParParams &);
    
    double estimate_error();
    virtual void update_data(double val);
    
protected:
    
    double f;
    double f2;
    
    
};

#endif	/* SIMPLEVAR_H */

