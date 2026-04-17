#include "criticalCells.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";

    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string promptRequiredLine(const std::string& prompt)
{
    while (true)
    {
        std::cout << prompt;

        std::string line;
        if (!std::getline(std::cin, line))
            throw std::runtime_error("Input stream closed.");

        line = trim(line);
        if (!line.empty())
            return line;

        std::cout << "  Please enter a value.\n";
    }
}

std::string promptFilename()
{
    while (true)
    {
        const std::string filename = promptRequiredLine("Input CSV filename: ");
        const std::filesystem::path filepath(filename);

        if (!std::filesystem::exists(filepath))
        {
            std::cout << "  File not found. Please enter an existing path.\n";
            continue;
        }

        if (!std::filesystem::is_regular_file(filepath))
        {
            std::cout << "  That path is not a regular file.\n";
            continue;
        }

        std::ifstream file_stream(filepath);
        if (!file_stream.good())
        {
            std::cout << "  Unable to open that file for reading.\n";
            continue;
        }

        return filename;
    }
}

size_t promptMaxDimension()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input max dimension: ");
        std::stringstream parser(line);
        long long value = 0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0)
        {
            std::cout << "  Please enter a positive integer.\n";
            continue;
        }

        return static_cast<size_t>(value);
    }
}

std::vector<double> promptEpsilonBreaks()
{
    while (true)
    {
        const std::string raw_line = promptRequiredLine("Input epsilon breaks (space or comma separated): ");
        std::string normalized = raw_line;
        std::replace(normalized.begin(), normalized.end(), ',', ' ');

        std::stringstream parser(normalized);
        std::vector<double> eps_breaks;
        std::string token;
        bool valid = true;

        while (parser >> token)
        {
            size_t parsed_length = 0;
            double value = 0.0;

            try
            {
                value = std::stod(token, &parsed_length);
            }
            catch (const std::exception&)
            {
                valid = false;
                break;
            }

            if (parsed_length != token.size() || value <= 0.0)
            {
                valid = false;
                break;
            }

            eps_breaks.push_back(value);
        }

        if (!valid || eps_breaks.empty())
        {
            std::cout << "  Please enter one or more positive numbers, for example: 1.0 1.4 1.7\n";
            continue;
        }

        bool strictly_increasing = true;
        for (size_t i = 1; i < eps_breaks.size(); ++i)
        {
            if (eps_breaks[i] <= eps_breaks[i - 1])
            {
                strictly_increasing = false;
                break;
            }
        }

        if (!strictly_increasing)
        {
            std::cout << "  Epsilon breaks must be strictly increasing.\n";
            continue;
        }

        return eps_breaks;
    }
}

int promptThreadNumber()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input thread number: ");
        std::stringstream parser(line);
        int value = 0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0)
        {
            std::cout << "  Please enter a positive integer.\n";
            continue;
        }

        return value;
    }
}

double promptPVCapScale()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input PV diameter cap scale: ");
        std::stringstream parser(line);
        double value = 0.0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0.0)
        {
            std::cout << "  Please enter a positive number.\n";
            continue;
        }

        return value;
    }
}

bool promptConfirmation(const std::string& filename, const size_t maxdim, const std::vector<double>& eps_breaks, const int threadnumber, const double pv_cap_scale)
{
    while (true)
    {
        std::cout << "\nRun morsePiecewisePH with:\n";
        std::cout << "  input file: " << filename << '\n';
        std::cout << "  max dimension: " << maxdim << '\n';
        std::cout << "  epsilon breaks:";
        for (const double eps : eps_breaks)
            std::cout << ' ' << eps;
        std::cout << '\n';
        std::cout << "  thread number: " << threadnumber << '\n';
        std::cout << "  pv cap scale: " << pv_cap_scale << '\n';

        const std::string answer = promptRequiredLine("Proceed? [y/n]: ");
        if (answer == "y" || answer == "Y" || answer == "yes" || answer == "YES" || answer == "Yes")
            return true;
        if (answer == "n" || answer == "N" || answer == "no" || answer == "NO" || answer == "No")
            return false;

        std::cout << "  Please answer y or n.\n";
    }
}
}

int main()
{
    try
    {
        std::cout << "morsePiecewisePH interactive runner\n\n";

        const std::string filename = promptFilename();
        const size_t maxdim = promptMaxDimension();
        const std::vector<double> eps_breaks = promptEpsilonBreaks();
        const int threadnumber = promptThreadNumber();
        const double pv_cap_scale = promptPVCapScale();

        if (!promptConfirmation(filename, maxdim, eps_breaks, threadnumber, pv_cap_scale))
        {
            std::cout << "Run cancelled.\n";
            return 0;
        }

        CritCells<VR, NormalDistMat> cc(filename);

        std::cout << "\nstarted running...\n";

        const auto st0 = std::chrono::high_resolution_clock::now();
        cc.morsePiecewisePH(maxdim, eps_breaks, threadnumber, pv_cap_scale);
        const auto st1 = std::chrono::high_resolution_clock::now();

        const auto pt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(st1 - st0);
        std::cout << "run time = " << pt_ms.count() << '\n';
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}

// int main()
// {   
//     std::string filename = "test_4dsphere_100.csv";
//     CritCells<VR, NormalDistMat> cc(filename);

//     std::cout<<"started running...\n";

//     auto st0 = std::chrono::high_resolution_clock::now();

//     // cc.morseQuotientAndExpand(3, 1.0, 1.7, 4);
//     // cc.morseVRPH(3, 1.4, 4);

//     std::vector<double> eps_breaks = {1.1, 1.4, 1.7};
//     cc.morsePiecewisePH(4, eps_breaks, 4, 0.8);

//     auto st1 = std::chrono::high_resolution_clock::now();
//     auto pt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(st1 - st0);
//     std::cout<<"run time = "<<pt_ms.count() <<'\n';

//     return 0;
// }

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
