#include "criticalCells.cpp"
#include <iostream>
#include <cstdlib> // for std::stoi

int main(int argc, char *argv[])
{
    switch (argc)
    {
    case 3:
        {
            CritCells<VR, NormallDistMat> cc(argv[1]);
            cc.run_Compute(std::stoi(argv[2]));
            break;
        }
    case 4:
        {
            CritCells<VR, NormallDistMat> cc(argv[1]);
            cc.run_Compute(std::stoi(argv[2]), std::stoi(argv[3]));
            break;
        }
    default:
        std::cerr << "Usage: " << argv[0] << " <input_filename> <maxDim> [<batch_size>]" << std::endl;
        return 1;
    }

    return 0;
}
