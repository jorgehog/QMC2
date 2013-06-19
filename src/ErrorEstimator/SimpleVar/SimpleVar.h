/* 
 * File:   SimpleVar.h
 * Author: jorgmeister
 *
 * Created on October 30, 2012, 6:13 PM
 */

#ifndef SIMPLEVAR_H
#define	SIMPLEVAR_H

#include "../ErrorEstimator.h"
struct ParParams;

/*! \brief Calculates the simple variance of the sampled values.
 */
class SimpleVar : public ErrorEstimator {
public:
    SimpleVar(ParParams &);
    
    double estimate_error();
    
    /*!
     * Overrides the default described in the superclass.
     * Does not store values in memory, but rather use sum variables.
     */
    void update_data(double val);
    
    void normalize(); //NOT NEEDED?
protected:
    
    double f; //!< sum variable used to calculate the mean
    double f2; //!< sum variable used to calulate the mean of squares.
    
    
};

#endif	/* SIMPLEVAR_H */

