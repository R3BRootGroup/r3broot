#pragma once
#include <Math/GenVector/Cartesian3D.h>
#include <Math/GenVector/Polar3D.h>
#include <Math/GenVector/PxPyPzE4D.h>
#include <TRandom3.h>
#include <cmath>
#include <fmt/format.h>

template <typename AngleDist, typename EnergyDist, typename PointDist>
class TrackGenerator
{
  public:
    TrackGenerator(AngleDist angle_dist, EnergyDist energy_dist, PointDist point_dist)
        : angle_dist_(angle_dist)
        , energy_dist_(energy_dist)
        , point_dist_(point_dist)
    {
    }

    // void set_detector_width(double width) { detector_width_ = width; }
    // void set_detector_length(double length) { detector_length_ = length; }
    // void set_detector_height(double height) { detector_height_ = height; }
    void set_detector_size(double detector_size) { detector_size_ = detector_size; }
    void set_particle_mass(double particle_mass) { particle_mass_ = particle_mass; }

  private:
    // double detector_width_{};
    // double detector_length_{};
    // double detector_height_{};
    double particle_mass_{};
    double detector_size_{ 50.0 };
    static constexpr double speed_of_light_ = 299792458.0;

    ROOT::Math::Cartesian3D<double> position_{};

    AngleDist angle_dist_;
    EnergyDist energy_dist_;
    PointDist point_dist_;
    TRandom3 rd_engine_;

    auto rd_num_gen_point() -> ROOT::Math::Cartesian3D<double>;
    auto rd_num_gen_angles() -> ROOT::Math::Polar3D<double>;
    auto rd_num_gen_energy() -> double;
    auto calculate_abs_momentum(double kinetic_energy) -> double;
    auto calculate_momentum_energy(double kinetic_energy, double Theta, double Phi) -> ROOT::Math::PxPyPzE4D<double>;

    void calculate_position();
};

template <typename AngleDist, typename EnergyDist, typename PointDist>
auto TrackGenerator<AngleDist, EnergyDist, PointDist>::rd_num_gen_angles() -> ROOT::Math::Polar3D<double>
{
    ROOT::Math::Polar3D<double> angles{};
    angles.SetPhi(rd_engine_.Uniform(0., M_PI));
    angles.SetTheta(angle_dist_(rd_engine_));
    return angles;
}

template <typename AngleDist, typename EnergyDist, typename PointDist>
auto TrackGenerator<AngleDist, EnergyDist, PointDist>::rd_num_gen_energy() -> double
{
    return energy_dist_(rd_engine_);
}

template <typename AngleDist, typename EnergyDist, typename PointDist>
auto TrackGenerator<AngleDist, EnergyDist, PointDist>::rd_num_gen_point() -> ROOT::Math::Cartesian3D<double>
{

    return point_dist_(rd_engine_);
}

template <typename AngleDist, typename EnergyDist, typename PointDist>
auto TrackGenerator<AngleDist, EnergyDist, PointDist>::calculate_abs_momentum(double kinetic_energy) -> double
{
    return std::sqrt(kinetic_energy * kinetic_energy +
                     2 * particle_mass_ * speed_of_light_ * speed_of_light_ * kinetic_energy);
}

template <typename AngleDist, typename EnergyDist, typename PointDist>
auto TrackGenerator<AngleDist, EnergyDist, PointDist>::calculate_momentum_energy(double kinetic_energy,
                                                                                 double Theta,
                                                                                 double Phi)
    -> ROOT::Math::PxPyPzE4D<double>
{
    ROOT::Math::PxPyPzE4D<double> momentum_energy{ 0, 0, 0, kinetic_energy };
    double abs_momentum = calculate_abs_momentum(kinetic_energy);
    momentum_energy.SetPx(abs_momentum * std::sin(Theta) * std::cos(Phi));
    momentum_energy.SetPy(abs_momentum * std::sin(Theta) * std::sin(Phi));
    momentum_energy.SetPz(abs_momentum * std::cos(Theta));
    return momentum_energy;
}

template <typename AngleDist, typename EnergyDist, typename PointDist>
void TrackGenerator<AngleDist, EnergyDist, PointDist>::calculate_position()
{
    auto point = rd_num_gen_point();
    auto angles = rd_num_gen_angles();
    auto energy = rd_num_gen_energy();
    auto momentum_energy = calculate_momentum_energy(energy, angles.Theta(), angles.Phi());

    if (momentum_energy.Pz() == 0)
    {

        position_.SetX(point.x() - std::sin(angles.Phi()) * detector_size_);
        position_.SetY(point.y() - std::cos(angles.Phi()) * detector_size_);
        position_.SetZ(point.z());
    }
    else
    {
        position_.SetX(point.x() - std::sin(angles.Theta()) * std::cos(angles.Phi()) * detector_size_);
        position_.SetY(point.y() - std::sin(angles.Theta()) * std::sin(angles.Phi()) * detector_size_);
        position_.SetZ(point.z() - std::cos(angles.Theta()) * detector_size_);
    }
}
