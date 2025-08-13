#include "criticalCells.hpp"
#include <iostream>
#include <cstdlib> // for std::stoi
#include <fstream>

#include <chrono>

int main()
{   
    std::string filename = "test_5d_60_double.csv";
    // std::string filename = "test_2dsquare_double_15.csv";
    CritCells<VR, NormalDistMat> cc(filename);

    // cc.runAlphaTest(filename);
    // auto critical_weight = cc.run_MorseMatchPersistence(3, 0, 17);

    // auto critical_weight = cc.run_MorseMatch(3, 0, 12, 1);

    // for(auto& dim_cw: critical_weight)
    // {
    //     for (auto& w: dim_cw) std::cout<<w<<"  ";
    //     std::cout<<'\n';./
    // }
    auto st0 = std::chrono::high_resolution_clock::now();
    // cc.runMorseTest(5, 12.0, 4);
    // cc.runAlphaTest(filename, 2.0, 4);
    cc.runMorseReductionTest(5, 2.66, 9, 4);
    auto st1 = std::chrono::high_resolution_clock::now();
    auto pt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(st1 - st0);
    std::cout<<"run time = "<<pt_ms.count() <<'\n';

    return 0;
}

// int main(int argc, char *argv[])
// {
//     std::ofstream logFile("logfile.txt", std::ios::app);
//     std::streambuf *originalClogBuffer = std::clog.rdbuf();
//     std::clog.rdbuf(logFile.rdbuf());

//     switch (argc)
//     {
//     case 3:
//     {
//         CritCells<VR, NormallDistMat> cc(argv[1]);
//         cc.run_Compute(std::stoi(argv[2]));
//         break;
//     }
//     case 4:
//     {
//         CritCells<VR, NormallDistMat> cc(argv[1]);
//         cc.run_Compute(std::stoi(argv[2]), std::stoi(argv[3]));
//         break;
//     }
//     default:
//         std::cerr << "Usage: " << argv[0] << " <input_filename> <maxDim> [<batch_size>]" << std::endl;
//         return 1;
//     }

//     std::clog.rdbuf(originalClogBuffer);
//     return 0;
// }
