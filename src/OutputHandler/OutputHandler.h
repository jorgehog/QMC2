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
    
    bool parallel;
    int my_rank;
    int num_procs;
    
    std::string filename;
    std::string path;

    std::ofstream file;

    QMC* qmc;

public:

    OutputHandler();
    OutputHandler(std::string filename,
            std::string path,
            bool parallel,
            int my_rank,
            int num_procs
            );

    virtual void dump() = 0;
    virtual void finalize();

    void set_qmc_ptr(QMC* qmc) {
        this->qmc = qmc;
    }

};


#endif	/* OUTPUTHANDLER_H */