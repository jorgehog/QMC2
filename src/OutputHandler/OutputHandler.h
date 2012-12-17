/* 
 * File:   OutputHandler.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef OUTPUTHANDLER_H
#define	OUTPUTHANDLER_H

class OutputHandler {
protected:
    bool is_vmc;
    bool is_dmc;

    bool is_ASGD;

    bool parallel;
    int node;
    int n_nodes;

    bool use_file;

    std::stringstream s;

    std::string filename;
    std::string path;

    std::ofstream file;

    QMC* qmc;
    DMC* dmc;
    VMC* vmc;
    Minimizer* min;
    ASGD* asgd;

    void init_file();

    virtual void post_pointer_init() {
    };

public:

    OutputHandler();
    OutputHandler(std::string filename,
            std::string path,
            bool parallel,
            int node,
            int n_nodes
            );

    virtual void dump() = 0;
    virtual void finalize();

    void set_qmc_ptr(QMC* qmc);
    void set_min_ptr(Minimizer* min);

};


#endif	/* OUTPUTHANDLER_H */