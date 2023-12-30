#include "criticalCells.cpp"

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_filename> <maxDim>" << std::endl;
        return 1;
    }
    CritCells<VR, NormallDistMat> cc(argv[1]);
    cc.run_Compute(std::stoi(argv[2]));
    return 0;
}