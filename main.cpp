#include "criticalCells.hpp"
#include <iostream>
#include <cstdlib> // for std::stoi
#include <fstream>

#include <chrono>

int main()
{   
    std::string filename = "test_3dsphere_100.csv";
    CritCells<VR, NormalDistMat> cc(filename);

    std::cout<<"started running...\n";

    auto st0 = std::chrono::high_resolution_clock::now();

    // cc.morseQuotientAndExpand(3, 1.3, 1.7, 4);
    // cc.morseVRPH(3, 1.7, 4);

    std::vector<double> eps_breaks = {0.9, 1.3, 1.7};
    cc.morsePiecewisePH(3, eps_breaks, 4);

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
