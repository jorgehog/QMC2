/* 
 * File:   BlockingData.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef BLOCKING_H
#define	BLOCKING_H

class BlockingData : public OutputHandler {
public:

    BlockingData(std::string filename = "blockdata_out",
            std::string path = "./",
            bool parallel = false,
            int my_rank = 0,
            int num_procs = 1
            );

    virtual void dump();

};


#endif	/* BLOCKING_H */