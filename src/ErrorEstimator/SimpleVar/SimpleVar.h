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
    SimpleVar();
    
    void update_data(double val){
        E += val;
        E2 += val*val;
        i += 1;
    }
    
    double estimate_error(){
        return (E2 - E*E/i)/i;
    }
    
protected:
    
    double E;
    double E2;
    
    
};

#endif	/* SIMPLEVAR_H */

