#include<libmfpic/VelocityHistogram.hpp>

namespace mfpic {

double VelocityHistogram::evaluate(double velocity) const
{
    const int nbins = static_cast<int>(bin_values.size());

    int k = static_cast<int>((velocity - vmin) / dv);

    k = std::max(0, std::min(k, nbins - 1));

    return bin_values[k];
}

void VelocityHistogram::writeToCSVFile(
    const std::string &base_filename,
    const int time_step,
    const double time) const
{
    std::string base = base_filename;
    const std::string ext = ".csv";
    if (base.size() >= ext.size() &&
        base.compare(base.size() - ext.size(), ext.size(), ext) == 0)
    {
        base.erase(base.size() - ext.size());
    }

    const std::string filename = base + ".csv";

    bool write_header = false;
    if (time_step == 0) write_header = true;

    std::ofstream csv_file(filename, std::ios::out | std::ios::app); 
    csv_file << std::setprecision(std::numeric_limits<double>::max_digits10);

    if (write_header) {
        csv_file << "# Time_Step Time Velocity f\n";
    }

    for (std::size_t k = 0; k < bin_values.size(); ++k)
    {
        const double velocity = 0.5 * (bin_edges[k] + bin_edges[k + 1]);

        csv_file
            << time_step << " "
            << time << " "
            << velocity << " "
            << bin_values[k] << "\n";
    }
}


} // namespace 