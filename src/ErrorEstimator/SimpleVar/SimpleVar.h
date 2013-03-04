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
    void update_data(double val);
    void normalize();
protected:
    
    double f;
    double f2;
    
    
};

#endif	/* SIMPLEVAR_H */

