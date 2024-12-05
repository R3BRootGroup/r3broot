#include "NeulandPointFilter.h"

void NeulandPointFilter::SetFilter(R3B::Neuland::BitSetParticle filtered_particles)
{
    filtered_particles_ = filtered_particles;
}
void NeulandPointFilter::SetFilter(R3B::Neuland::BitSetParticle filtered_particles, double minimum_allowed_energy)
{
    filtered_particles_ = filtered_particles;
    minimum_allowed_energy_ = minimum_allowed_energy;
}

auto NeulandPointFilter::CheckFiltered(const R3BNeulandPoint& neuland_point) -> bool
{
    return (
        R3B::Neuland::CheckCriteria(R3B::Neuland::PidToBitSetParticle(neuland_point.GetPID()), filtered_particles_) or
        (neuland_point.GetEnergyLoss() < minimum_allowed_energy_));
}
