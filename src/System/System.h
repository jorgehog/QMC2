/* 
 * File:   System.h
 * Author: jorgehog
 *
 * Created on 30. mars 2012, 16:49
 */

#ifndef SYSTEM_H
#define	SYSTEM_H

class System {
protected:
    int n_p;
    int n2;
    int dim;
    double a_sym;
    double a_asym;

    std::vector<Potential*> potentials;
    Orbitals *orbital;

public:
    System(int n_p, int dim, Orbitals* orbital);
    System();
    
    void add_potential(Potential* pot);

    double get_potential_energy(const Walker* walker);
    
    virtual void update_walker(Walker* walker_pre, const Walker* walker_post, int particle) const = 0;
    
    virtual void calc_for_newpos(const Walker* walker_old, Walker* walker_new, int particle) = 0;

    virtual double get_spatial_ratio(const Walker* walker_pre, const Walker* walker_post, int particle) const = 0;

    virtual double get_spatial_wf(const Walker* walker) = 0;
    virtual void get_spatial_grad(Walker* walker, int particle) const = 0;
    virtual double get_spatial_lapl_sum(const Walker* walker) const = 0;
    virtual void initialize(Walker* walker) = 0;

    virtual void copy_walker(const Walker* parent, Walker* child) const = 0;
    virtual void reset_walker_ISCF(const Walker* walker_pre, Walker* walker_post, int particle) const = 0;
    
    Orbitals* get_orbital_ptr() {
        return orbital;
    }
};

#endif	/* SYSTEM_H */
