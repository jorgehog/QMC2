/* 
 * File:   Distribution.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef DISTRIBUTION_H
#define	DISTRIBUTION_H

class Distribution : public OutputHandler {
public:
    Distribution(ParParams &, std::string path, std::string name);

    void dump();
    void finalize();

private:
    
    std::string name;
    void generate_distribution();

};


#endif	/* DISTRIBUTION_H */