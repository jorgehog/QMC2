/* 
 * File:   None.h
 * Author: jorgmeister
 *
 * Created on October 29, 2012, 3:59 PM
 */

#ifndef NONE_H
#define	NONE_H

class None : public ErrorEstimator {
public:
    None();

    double estimate_error() {
        return 0;
    }

    virtual void update_data(double val) {

    }

protected:

};

#endif	/* NONE_H */

