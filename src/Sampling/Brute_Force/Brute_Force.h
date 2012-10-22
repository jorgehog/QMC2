/* 
 * File:   Brute_Force.h
 * Author: jorgmeister
 *
 * Created on October 12, 2012, 2:43 PM
 */

#ifndef BRUTE_FORCE_H
#define	BRUTE_FORCE_H

class Brute_Force : public Sampling {
public:

    Brute_Force(GeneralParams &);

    void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const {
        
    }

    virtual void get_necessities(Walker* walker){
        
    }
    
    virtual void update_necessities(const Walker* walker_pre, Walker* walker_post, int particle) const {
        
    }

    virtual void calculate_energy_necessities(Walker* walker) const {
        qmc->get_gradients(walker);
    }

    virtual void copy_walker(const Walker* parent, Walker* child) const {
        
    }
    
    virtual void reset_walker(const Walker* walker_pre, Walker* walker_post, int particle) const {
        
    }

};

#endif	/* BRUTE_FORCE_H */

