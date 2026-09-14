#include<libmfpic/VelocityHistogram.hpp>

namespace mfpic {

void VelocityHistogram::buildVelocityHistogram(
    const ParticleContainer& particles,
    mfem::Mesh& mesh,
    int nbins)
{
    double vmin = std::numeric_limits<double>::max();
    double vmax = -std::numeric_limits<double>::max();

    for (const Particle& particle : particles) {
        if (!particle.is_alive) continue;

        double vx = particle.velocity(0);

        vmin = std::min(vmin, vx);
        vmax = std::max(vmax, vx);
    }

    vmin_ = vmin;
    vmax_ = vmax;
    dv_   = (vmax - vmin) / nbins;

    bin_values_.assign(nbins, 0.0);
    bin_edges_.resize(nbins + 1);

    for (int k = 0; k <= nbins; ++k)
        bin_edges_[k] = vmin + k * dv_;

    for (const Particle& particle : particles) {
        if (!particle.is_alive) continue;

        double vx = particle.velocity(0);

        int k = static_cast<int>((vx - vmin) / dv_);

        k = std::max(0, std::min(k, nbins - 1));

        bin_values_[k] += particle.weight;
    }

    double domain_volume = 0.0;
    for (int e = 0; e < mesh.GetNE(); ++e)
    {
        domain_volume += mesh.GetElementVolume(e);
    }

    //normalize to number density
    const double inverse_volume_dv= 1.0 / (domain_volume * dv_);
    for (auto& value : bin_values_)
        value *= inverse_volume_dv;

    //normalize to one
    // for (auto& value : bin_values_)
    //     value /= (total_weight * dv_);


}

double VelocityHistogram::evaluate(double velocity) const
{
    const int nbins = static_cast<int>(bin_values_.size());

    int k = static_cast<int>((velocity - vmin_) / dv_);

    k = std::max(0, std::min(k, nbins - 1));

    return bin_values_[k];
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

    for (std::size_t k = 0; k < bin_values_.size(); ++k)
    {
        const double velocity = 0.5 * (bin_edges_[k] + bin_edges_[k + 1]);

        csv_file
            << time_step << " "
            << time << " "
            << velocity << " "
            << bin_values_[k] << "\n";
    }
}


} // namespace 