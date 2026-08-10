#pragma once
#include <libmfpic/ParticleContainer.hpp>
#include <mfem/mfem.hpp>

namespace mfpic {

class VelocityHistogram
{
public:

    void buildVelocityHistogram(
        const ParticleContainer& particles,
        mfem::Mesh& mesh,
        int nbins
    );

    double evaluate(double velocity) const;

    void writeToCSVFile(
        const std::string &base_filename,
        int time_step,
        double time) const;
    
    double getVMin() const { return vmin_; }
    double getVMax() const { return vmax_; }
    double getDV()   const { return dv_;   }

    const std::vector<double>& getBinValues() const { return bin_values_; }
    const std::vector<double>& getBinEdges()  const { return bin_edges_;  }

    std::size_t getNumBins() const { return bin_values_.size(); }
    std::size_t getNumEdges() const { return bin_edges_.size(); }
    bool empty() const { return bin_values_.empty(); }

private:

    double vmin_ = 0.0;
    double vmax_ = 0.0;
    double dv_ = 0.0;

    std::vector<double> bin_values_;
    std::vector<double> bin_edges_;
};

} //namespace 