#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <exception>

#include "nwqec/gridsynth/gridsynth.hpp"
#include "nwqec/core/constants.hpp"

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <epsilon> <theta_input.txt> <output.csv>\n\n"
        << "theta_input.txt should contain one theta per line.\n";
}

struct GateCounts {
    int t = 0;
    int h = 0;
    int s = 0;
    int w = 0;
};

static GateCounts count_gates(const std::string& gates) {
    GateCounts c;
    for (char g : gates) {
        if (g == 'T') ++c.t;
        else if (g == 'H') ++c.h;
        else if (g == 'S') ++c.s;
        else if (g == 'W') ++c.w;
    }
    return c;
}

int main(int argc, char* argv[]) {
    // Expect: prog epsilon theta_input.txt output.csv
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string epsilon = argv[1];
    const std::string theta_file = argv[2];
    const std::string out_name = argv[3];

    // Open theta input file
    std::ifstream ifs(theta_file);
    if (!ifs.is_open()) {
        std::cerr << "Error: cannot open theta input file: " << theta_file << "\n";
        return 2;
    }

    // Open output file
    std::ofstream ofs(out_name, std::ios::out | std::ios::trunc);
    if (!ofs.is_open()) {
        std::cerr << "Error: cannot open output file: " << out_name << "\n";
        return 3;
    }

    // Write CSV header
    ofs << "theta,t_count,h_count,s_count,w_count\n";

    std::string theta;
    while (std::getline(ifs, theta)) {

        // Skip empty lines
        if (theta.empty())
            continue;

        try {
            const auto gates = gridsynth::gridsynth_gates(
                theta,
                epsilon,
                NWQEC::DEFAULT_DIOPHANTINE_TIMEOUT_MS,
                NWQEC::DEFAULT_FACTORING_TIMEOUT_MS,
                false,
                false
            );

            const GateCounts c = count_gates(gates);

            ofs << theta << ","
                << c.t << ","
                << c.h << ","
                << c.s << ","
                << c.w << "\n";

        } catch (const std::exception& e) {
            std::cerr << "Warning: failed for theta=" << theta
                      << " : " << e.what() << "\n";
            ofs << theta << ",-1,-1,-1,-1\n";
        } catch (...) {
            std::cerr << "Warning: failed for theta=" << theta
                      << " : unknown error\n";
            ofs << theta << ",-1,-1,-1,-1\n";
        }
    }

    ifs.close();
    ofs.close();

    std::cout << "Results written to: " << out_name << "\n";
    return 0;
}
