#ifndef GRAPHCMD_HELPER_H
#define GRAPHCMD_HELPER_H

#include <iostream>
#include <fstream>

#include "../common/cmd_helpers.h"
#include "binopt.h"

namespace cmd
{
    void saveAssignments(const std::string& filename ,
                         int n_sites,
                         BinaryOptimization<>& bin_opt)
    {
        using namespace std;

        ofstream out_file;
        out_file.open(filename, ios::trunc);

        /*
        for(auto row = 0; row < n_sites; row++) {
            out_file << row
                     << "\t"
                     << bin_opt.whichLabel(row + 2)
                     << endl;
        }
        */
        out_file.close();
    }
}

#endif // GRAPHCMD_HELPER_H