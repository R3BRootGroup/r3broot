#include "generator.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/Polar3D.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include <algorithm>
#include <cmath>
#include <fmt/format.h>
#include <TRandom.h>
template <typename ADist, typename EDist>
void TrackGenerator<ADist, EDist>::set_detector_width(double width)
{
    detector_width = width;
}
template <typename ADist, typename EDist>
void TrackGenerator<ADist, EDist>::set_detector_length(double length) { detector_length = length; }
template <typename ADist, typename EDist>
void TrackGenerator<ADist, EDist>::set_detector_height(double height) { detector_height = height; }

template <typename ADist, typename EDist>
void TrackGenerator<ADist, EDist>::calculate_postion()
{
    auto point = rd_num_gen_point();
    auto angles = rd_num_gen_angles();
    auto energy = rd_num_gen_energy();
    auto momentum_energy = calculate_momentum_energy(energy, angles.Theta(), angles.Phi());
    if (momentum_energy.Pz() == 0)
    {
        double radius = std::sqrt(detector_width * detector_width + detector_length * detector_length);

        position.SetX(point.x() - std::sin(angles.Phi()) * radius);

        position.SetY(point.y() - std::cos(angles.Phi()) * radius);

        position.SetZ(point.z());
    }
    else
    {
        double diagonal = std::sqrt(detector_width * detector_width + detector_length * detector_length +
                                    detector_height * detector_height);

        position.SetX(point.x() - std::sin(angles.Theta()) * std::cos(angles.Phi()) * diagonal);

        position.SetY(point.y() - std::sin(angles.Theta()) * std::sin(angles.Phi()) * diagonal);

        position.SetZ(point.z() - std::cos(angles.Theta()) * diagonal);
    }
}

template <typename ADist, typename EDist>
double TrackGenerator<ADist, EDist>::calculate_abs_momentum(double kinetic_energy)
{

    return std::sqrt(kinetic_energy * kinetic_energy +
                     2 * particle_mass * speed_of_light * speed_of_light * kinetic_energy);
}

template <typename ADist, typename EDist>
ROOT::Math::PxPyPzE4D<double> TrackGenerator<ADist, EDist>::calculate_momentum_energy(double kinetic_energy, double Theta, double Phi)
{
    ROOT::Math::PxPyPzE4D<double> momentum_energy{ 0, 0, 0, kinetic_energy };
    double abs_momentum = calculate_abs_momentum(kinetic_energy);
    momentum_energy.SetPx(abs_momentum * std::sin(Theta) * std::cos(Phi));
    momentum_energy.SetPy(abs_momentum * std::sin(Theta) * std::sin(Phi));
    momentum_energy.SetPz(abs_momentum * std::cos(Theta));
    return momentum_energy;
}

template <typename ADist, typename EDist>
ROOT::Math::Cartesian3D<double> TrackGenerator<ADist, EDist>::rd_num_gen_point()
{
    ROOT::Math::Cartesian3D<double> point{};
    
    point.SetX();
    // auto dist1 = std::uniform_real_distribution<>{ 0., detector_length };
    // auto dist2 = std::uniform_real_distribution<>{ 0., detector_width };
    // auto dist3 = std::uniform_real_distribution<>{ 0., detector_height };
    //
    //
    // point.SetX(dist1(gen));
    // point.SetY(dist2(gen));
    // point.SetZ(dist3(gen));
 return point;
} };
    //
    //
    // point.SetX(dist1(gen));
    // point.SetY(dist2(gen));
    // point.SetZ(dist3(gen));
 return point;
}
