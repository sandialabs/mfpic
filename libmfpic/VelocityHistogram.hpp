#pragma once
#include <mfem/mfem.hpp>

namespace mfpic {

class VelocityHistogram
{
public:
    double evaluate(double velocity) const;

    void writeToCSVFile(
        const std::string &base_filename,
        int time_step,
        double time) const;

    double vmin = 0.0;
    double vmax = 0.0;
    double dv = 0.0;

    std::vector<double> bin_values;
    std::vector<double> bin_edges;
};

} //namespace 