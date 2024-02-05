/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BNEULANDHITMON_H
#define R3BNEULANDHITMON_H

/** Neuland Digi Monitor
 * @author Jan Mayer
 * @since  10.06.2015
 *
 *  Input:  Digis. Should work both on Digis from Monte Carlo simulations
 *          as well as experimental data on digi level.
 *  Output: Vairous diagrams, 3D scatterplot for each event if option 3DTRACK is given.
 */

#include "FairTask.h"
#include "R3BDataMonitor.h"
#include "R3BIOConnector.h"
#include "R3BNeulandHit.h"
#include "TCAConnector.h"
#include <string>

class TH1D;
class TH2D;
class TH3D;
class TH1I;

class R3BNeulandHitMon : public FairTask
{
  public:
    explicit R3BNeulandHitMon(const Option_t* option = "");

    void Exec(Option_t* /*option*/) override;

    void SetDistanceToTarget(double distance) { distance_to_target_ = distance; }

  protected:
    auto Init() -> InitStatus override;
    void Finish() override { data_monitor_.save_to_sink(); }

  private:
    std::string output_{ "NeulandHitMon" };

    R3B::InputVectorConnector<R3BNeulandHit> neuland_hits_{ "NeulandHits" };

    double distance_to_target_ = 0.;

    bool is_3d_track_enabled_ = false;

    R3B::DataMonitor data_monitor_;

    TH3D* hist_3_ = nullptr;
    TH1D* hist_time_ = nullptr;
    TH1D* hist_time_adj_ = nullptr;
    TH1I* hist_mult_ = nullptr;
    TH1D* hist_depth_ = nullptr;
    TH1D* hist_foremost_energy_ = nullptr;
    TH1D* hist_sternmost_energy_ = nullptr;
    TH2D* hist_depth_vs_foremost_energy_ = nullptr;
    TH2D* hist_depth_vs_sternmost_energy_ = nullptr;
    TH1D* hist_energy_tot_ = nullptr;
    TH1D* hist_energy_ = nullptr;
    TH1D* hist_x_ = nullptr;
    TH1D* hist_y_ = nullptr;
    TH1D* hT = nullptr;
    TH1D* hTNeigh = nullptr;
    TH2D* hist_depth_vs_energy_tot_ = nullptr;
    TH2D* hist_pos_vs_energy_ = nullptr;
    TH2D* hdeltaEE = nullptr;
    TH1D* hist_beta_ = nullptr;

    ClassDefOverride(R3BNeulandHitMon, 0); // NOLINT
};

#endif // R3BNEULANDHITMON_H
