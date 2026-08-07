#include<libmfpic/VelocityHistogram.hpp>

namespace mfpic {

double VelocityHistogram::evaluate(double velocity) const
{
    const int nbins = static_cast<int>(bin_values.size());

    int k = static_cast<int>((velocity - vmin) / dv);

    k = std::max(0, std::min(k, nbins - 1));

    return bin_values[k];
}

} // namespace 