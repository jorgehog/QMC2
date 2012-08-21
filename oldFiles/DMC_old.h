/* 
 * File:   DMC.h
 * Author: jorgehog
 *
 * Created on 16. februar 2012, 16:18
 */


#ifndef DMC_H
#define	DMC_H



#include "Walker.h"
#include "lib.h"
#include "System.h"

class Walker;
class System;

class DMC {
protected:
    int n_c;
    int n_w;
    int dim;
    int n_p;
    int K; //how many times the walkers multiply their population;
    int step;
    bool importance;
    double D;
    double std;
    double dt;
    double E_T; //Trial energy
    double E;
    double GB;
 
public:
    int n_w_now;
    int n_w_orig;
    long ranseed;
    Walker **Angry_mob;
    Walker *current_walker;
    System* system;
    int* population;

    DMC();
    ~DMC();
    virtual void do_DMC() = 0;
  
    void increase_walker_space();
    void calc_gs_statistics();
    void bury_the_dead();
    void Evolve_walker(int k);
    void copy_walker(int parent, int child);
    bool singular(int k);
    virtual void update_pos(int k) = 0;
    int get_dim() const;
    int get_N_P() const;
    virtual bool is_importance() const = 0;
    
    virtual void initialize_walkers();
    virtual double get_new_pos() = 0;
};

//////////////
//BRUTE FORCE
/////////////

class DMC_BF : public DMC {
protected:
    double **tmp_r;
public:

    DMC_BF();
    ~DMC_BF();
    DMC_BF(double D, long randseed, int N_C, int N_W, System* system, double dt, double E_T);
    virtual void do_DMC();
    virtual bool is_importance() const;
    virtual void update_pos(int k);
    virtual double get_new_pos();

};


//////////////
//IMPORTANCE SAMPLING 
/////////////


class DMC_Importance : public DMC{
 public:
  DMC_Importance();
};






#endif	/* DMC_H */

