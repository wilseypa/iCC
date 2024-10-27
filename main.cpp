#include "criticalCells.cpp"
#include <iostream>
#include <cstdlib> // for std::stoi
#include <fstream>

int main()
{   std::string filename = "testdata3d.csv";
    CritCells<VR, NormallDistMat> cc(filename);

    auto critical_weight = cc.run_MorseMatch(3, 0, 8);

    // for(auto& dim_cw: critical_weight)
    // {
    //     for (auto& w: dim_cw) std::cout<<w<<"  ";
    //     std::cout<<'\n';
    // }

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
