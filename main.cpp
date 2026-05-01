#include "criticalCells.hpp"

#include <CLI/CLI.hpp>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
constexpr const char* PIECEWISE_PH_TOOL = "piecewise";
constexpr const char* MORSE_PH_TOOL = "ph";

class HelpSpacingFormatter : public CLI::Formatter
{
public:
    std::string make_option(const CLI::Option* option, bool is_positional) const override
    {
        std::string text = CLI::Formatter::make_option(option, is_positional);
        const auto& long_names = option->get_lnames();

        if (!is_positional &&
            std::find(long_names.begin(), long_names.end(), "help") != long_names.end())
        {
            text += '\n';
        }

        return text;
    }
};

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";

    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string normalizeTool(std::string tool)
{
    tool = trim(tool);
    std::transform(tool.begin(), tool.end(), tool.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (tool == "1" || tool == "piecewise" || tool == "morsepiecewiseph" || tool == "morse-piecewise-ph")
        return PIECEWISE_PH_TOOL;

    if (tool == "2" || tool == "ph" || tool == "morseph" || tool == "morse-ph")
        return MORSE_PH_TOOL;

    throw std::invalid_argument("Unknown tool '" + tool + "'. Use 'piecewise' or 'ph'.");
}

bool isPiecewisePHTool(const std::string& tool)
{
    return normalizeTool(tool) == PIECEWISE_PH_TOOL;
}

const char* displayToolName(const std::string& tool)
{
    return isPiecewisePHTool(tool) ? "morsePiecewisePH" : "morsePH";
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

std::string promptTool()
{
    while (true)
    {
        std::cout << "Select tool:\n";
        std::cout << "  1) morsePiecewisePH\n";
        std::cout << "  2) morsePH\n";

        const std::string answer = promptRequiredLine("Input tool [1/2]: ");
        try
        {
            return normalizeTool(answer);
        }
        catch (const std::invalid_argument&)
        {
            std::cout << "  Please enter 1 for morsePiecewisePH or 2 for morsePH.\n";
        }
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

double promptMaxEpsilon()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input max epsilon: ");
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

void printRunConfiguration(const std::string& tool,
                           const std::string& filename,
                           const size_t maxdim,
                           const int threadnumber,
                           const std::vector<double>& eps_breaks,
                           const double maxeps,
                           const double pv_cap_scale)
{
    const std::string normalized_tool = normalizeTool(tool);

    std::cout << "\nRun " << displayToolName(normalized_tool) << " with:\n";
    std::cout << "  input file: " << filename << '\n';
    std::cout << "  max dimension: " << maxdim << '\n';
    std::cout << "  thread number: " << threadnumber << '\n';

    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        std::cout << "  epsilon breaks:";
        for (const double eps : eps_breaks)
            std::cout << ' ' << eps;
        std::cout << '\n';
        std::cout << "  pv cap scale: " << pv_cap_scale << '\n';
    }
    else
    {
        std::cout << "  max epsilon: " << maxeps << '\n';
    }
}

bool promptConfirmation(const std::string& tool,
                        const std::string& filename,
                        const size_t maxdim,
                        const int threadnumber,
                        const std::vector<double>& eps_breaks,
                        const double maxeps,
                        const double pv_cap_scale)
{
    while (true)
    {
        printRunConfiguration(tool, filename, maxdim, threadnumber, eps_breaks, maxeps, pv_cap_scale);

        const std::string answer = promptRequiredLine("Proceed? [y/n]: ");
        if (answer == "y" || answer == "Y" || answer == "yes" || answer == "YES" || answer == "Yes")
            return true;
        if (answer == "n" || answer == "N" || answer == "no" || answer == "NO" || answer == "No")
            return false;

        std::cout << "  Please answer y or n.\n";
    }
}

void validateEpsilonBreaks(const std::vector<double>& eps_breaks)
{
    if (eps_breaks.empty())
        throw std::invalid_argument("--eps-breaks requires one or more positive numbers when --tool piecewise is selected.");

    for (const double eps : eps_breaks)
    {
        if (eps <= 0.0)
            throw std::invalid_argument("--eps-breaks values must be positive.");
    }

    for (size_t i = 1; i < eps_breaks.size(); ++i)
    {
        if (eps_breaks[i] <= eps_breaks[i - 1])
            throw std::invalid_argument("--eps-breaks values must be strictly increasing.");
    }
}

void validateRunOptions(const std::string& tool,
                        const std::string& filename,
                        const size_t maxdim,
                        const int threadnumber,
                        const std::vector<double>& eps_breaks,
                        const double maxeps,
                        const double pv_cap_scale)
{
    const std::string normalized_tool = normalizeTool(tool);

    if (filename.empty())
        throw std::invalid_argument("Input filename is required.");

    const std::filesystem::path filepath(filename);
    if (!std::filesystem::exists(filepath))
        throw std::invalid_argument("Input file does not exist: " + filename);
    if (!std::filesystem::is_regular_file(filepath))
        throw std::invalid_argument("Input path is not a regular file: " + filename);

    if (maxdim == 0)
        throw std::invalid_argument("Maximum dimension must be a positive integer.");
    if (threadnumber <= 0)
        throw std::invalid_argument("Thread number must be a positive integer.");

    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        validateEpsilonBreaks(eps_breaks);
        if (pv_cap_scale <= 0.0)
            throw std::invalid_argument("--pv-cap-scale is required and must be positive when --tool piecewise is selected.");
    }
    else
    {
        if (maxeps <= 0.0)
            throw std::invalid_argument("--max-eps is required and must be positive when --tool ph is selected.");
    }
}

int runTool(const std::string& tool,
            const std::string& filename,
            const size_t maxdim,
            const int threadnumber,
            const std::vector<double>& eps_breaks,
            const double maxeps,
            const double pv_cap_scale,
            const bool print_configuration,
            const bool verbose)
{
    const std::string normalized_tool = normalizeTool(tool);

    validateRunOptions(normalized_tool, filename, maxdim, threadnumber, eps_breaks, maxeps, pv_cap_scale);

    if (print_configuration && verbose)
        printRunConfiguration(normalized_tool, filename, maxdim, threadnumber, eps_breaks, maxeps, pv_cap_scale);

    CritCells<VR, NormalDistMat> icc(filename);

    if (verbose)
        std::cout << "\nstarted running...\n";

    const auto st0 = std::chrono::high_resolution_clock::now();
    if (normalized_tool == PIECEWISE_PH_TOOL)
        icc.morsePiecewisePH(maxdim, eps_breaks, threadnumber, pv_cap_scale, verbose);
    else
        icc.morseVRPH(maxdim, maxeps, threadnumber);
    const auto st1 = std::chrono::high_resolution_clock::now();

    const auto pt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(st1 - st0);
    if (verbose)
        std::cout << "run time = " << pt_ms.count() << '\n';

    return 0;
}

int runInteractive()
{
    std::cout << "Morse PH interactive runner\n\n";

    const std::string filename = promptFilename();
    const std::string tool = promptTool();
    const size_t maxdim = promptMaxDimension();
    const int threadnumber = promptThreadNumber();

    std::vector<double> eps_breaks;
    double maxeps = 0.0;
    double pv_cap_scale = 0.0;

    if (isPiecewisePHTool(tool))
    {
        eps_breaks = promptEpsilonBreaks();
        pv_cap_scale = promptPVCapScale();
    }
    else
    {
        maxeps = promptMaxEpsilon();
    }

    if (!promptConfirmation(tool, filename, maxdim, threadnumber, eps_breaks, maxeps, pv_cap_scale))
    {
        std::cout << "Run cancelled.\n";
        return 0;
    }

    return runTool(tool, filename, maxdim, threadnumber, eps_breaks, maxeps, pv_cap_scale, false, true);
}

int runCommandLine(int argc, char** argv)
{
    CLI::App app{"Morse PH command-line runner"};
    app.formatter(std::make_shared<HelpSpacingFormatter>());

    std::string tool;
    std::string filename;
    size_t maxdim = 0;
    int threadnumber = 1;
    std::vector<double> eps_breaks;
    double maxeps = 0.0;
    double pv_cap_scale = 0.0;
    bool verbose = false;

    app.add_option("-t,--tool", tool, "Tool to run: piecewise or ph")
        ->required();

    app.add_option("-f,--file-name", filename, "Input CSV filename")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("-d,--max-dim", maxdim, "Maximum simplex dimension")
        ->required()
        ->check(CLI::PositiveNumber);

    app.add_option("-n,--threads", threadnumber, "Number of worker threads")
        ->default_val(1)
        ->check(CLI::PositiveNumber);

    app.add_flag("-v,--verbose", verbose, "Show diagnostic output");

    // For a std::vector option, CLI11 already accepts one or more values when
    // the option is present. Conditional requirement and strictly-increasing
    // validation are handled after parsing in validateRunOptions().
    app.add_option("--eps-breaks", eps_breaks,
                   "Required when --tool piecewise. Strictly increasing epsilon breaks. Accepts spaces or commas, e.g. '--eps-breaks 1.0 1.4 1.7' or '--eps-breaks=1.0,1.4,1.7'")
        ->delimiter(',')
        ->check(CLI::PositiveNumber);

    app.add_option("--pv-cap-scale", pv_cap_scale,
                   "Required when --tool piecewise. PV diameter cap scale")
        ->check(CLI::PositiveNumber);

    app.add_option("--max-eps", maxeps,
                   "Required when --tool ph. Maximum epsilon")
        ->check(CLI::PositiveNumber);

    try
    {
        app.parse(argc, argv);
        return runTool(tool, filename, maxdim, threadnumber, eps_breaks, maxeps, pv_cap_scale, true, verbose);
    }
    catch (const CLI::ParseError& err)
    {
        return app.exit(err);
    }
}
}

int main(int argc, char** argv)
{
    try
    {
        // Backward-compatible behavior: running with no arguments keeps the old interactive front end.
        if (argc == 1)
            return runInteractive();

        return runCommandLine(argc, argv);
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}

// int main()
// {   
//     std::string filename = "test_4dsphere_100.csv";
//     CritCells<VR, NormalDistMat> icc(filename);

//     std::cout<<"started running...\n";

//     auto st0 = std::chrono::high_resolution_clock::now();

//     // icc.morseQuotientAndExpand(3, 1.0, 1.7, 4);
//     // icc.morseVRPH(4, 1.7, 4);

//     std::vector<double> eps_breaks = {1.4, 1.7};
//     icc.morsePiecewisePH(4, eps_breaks, 4, 0.8);

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
//         CritCells<VR, NormallDistMat> icc(argv[1]);
//         icc.run_Compute(std::stoi(argv[2]));
//         break;
//     }
//     case 4:
//     {
//         CritCells<VR, NormallDistMat> icc(argv[1]);
//         icc.run_Compute(std::stoi(argv[2]), std::stoi(argv[3]));
//         break;
//     }
//     default:
//         std::cerr << "Usage: " << argv[0] << " <input_filename> <maxDim> [<batch_size>]" << std::endl;
//         return 1;
//     }

//     std::clog.rdbuf(originalClogBuffer);
//     return 0;
// }
