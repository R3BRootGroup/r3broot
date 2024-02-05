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

#pragma once

#include "FairTask.h"
#include "Filterable.h"
#include "NeulandPointFilter.h"
#include "R3BDataMonitor.h"
#include "R3BDigitizingEngine.h"
#include "R3BDigitizingPaddleNeuland.h"
#include "R3BDigitizingTacQuila.h"
#include "R3BDigitizingTamex.h"
#include "R3BIOConnector.h"
#include "R3BNeulandGeoPar.h"
#include "R3BNeulandHit.h"
#include "R3BNeulandHitPar.h"
#include "R3BNeulandPoint.h"
#include "TCAConnector.h"
#include <TClonesArray.h>
#include <TH1.h>

class TGeoNode;
class TH1F;
class TH2F;

/**
 * NeuLAND digitizing finder task
 * @author Jan Mayer
 *
 * For each event, get the R3BNeulandPoints and apply the detector response to create R3BNeulandHits.
 *   Input:  Branch NeulandPoints = TClonesArray("R3BNeulandPoint")
 *           Stored Neuland Geometry Parameter NeulandGeoPar
 *   Output: Branch NeulandHits  = TClonesArray("R3BNeulandDigi")
 *   Additional output: Some control histograms
 *
 */
namespace Digitizing = R3B::Digitizing;

class R3BNeulandDigitizer : public FairTask
{
  public:
    enum class Options
    {
        neulandTamex,
        neulandTacquila
    };
    using NeulandPaddle = Digitizing::Neuland::NeulandPaddle;
    using TacquilaChannel = Digitizing::Neuland::TacQuila::Channel;
    using TamexChannel = Digitizing::Neuland::Tamex::Channel;
    template <typename Type>
    using UseChannel = Digitizing::UseChannel<Type>;
    template <typename Type>
    using UsePaddle = Digitizing::UsePaddle<Type>;

    R3BNeulandDigitizer();
    explicit R3BNeulandDigitizer(std::unique_ptr<Digitizing::DigitizingEngineInterface> engine);

  protected:
    auto Init() -> InitStatus override;
    void Finish() override { data_monitor_.save_to_sink(); }
    void SetParContainers() override;

  public:
    void Exec(Option_t* /*option*/) override;
    void SetEngine(std::unique_ptr<Digitizing::DigitizingEngineInterface> engine);
    void AddFilter(const Filterable<R3BNeulandHit&>::Filter& filter) { neuland_hit_filters_.Add(filter); }
    void SetNeulandPointFilter(R3B::Neuland::BitSetParticle particle);
    void SetNeulandPointFilter(R3B::Neuland::BitSetParticle particle, double minimum_allowed_energy_gev);

  private:
    R3B::InputVectorConnector<R3BNeulandPoint> neuland_points_{ "NeulandPoints" };
    R3B::OutputVectorConnector<R3BNeulandHit> neuland_hits_{ "NeulandHits" };

    std::unique_ptr<Digitizing::DigitizingEngineInterface> digitizing_engine_; // owning

    Filterable<R3BNeulandHit&> neuland_hit_filters_;

    R3BNeulandGeoPar* neuland_geo_par_ = nullptr; // non-owning
    NeulandPointFilter neuland_point_filter_;

    R3B::DataMonitor data_monitor_;
    TH1I* mult_one_ = nullptr;
    TH1I* mult_two_ = nullptr;
    TH1F* rl_time_to_trig_ = nullptr;

  public:
    template <typename... Args>
    explicit R3BNeulandDigitizer(Options option, Args&&... args)
        : R3BNeulandDigitizer()
    {
        switch (option)
        {
            case Options::neulandTamex:
                digitizing_engine_ = Digitizing::CreateEngine(UsePaddle<NeulandPaddle>(),
                                                              UseChannel<TamexChannel>(std::forward<Args>(args)...));
                break;
            case Options::neulandTacquila:
                digitizing_engine_ = Digitizing::CreateEngine(UsePaddle<NeulandPaddle>(),
                                                              UseChannel<TacquilaChannel>(std::forward<Args>(args)...));
                break;
        }
    }
    ClassDefOverride(R3BNeulandDigitizer, 1) // NOLINT
};
