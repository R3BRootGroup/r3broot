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

#ifndef R3BNEULAND_H
#define R3BNEULAND_H

#include "R3BDetector.h"
#include "R3BIOConnector.h"
#include "R3BNeulandPoint.h"
#include "TLorentzVector.h"
#include <memory>
#include <string>

class FairVolume;
class TClonesArray;
class R3BNeulandGeoPar;

/**
 * NeuLAND detector simulation class
 * @author Jan Mayer
 * @since 12.01.2016
 * For each simulated event, TClonesArrays are filled:
 * - NeulandPoints (R3BNeulandPoint), each representing energy deposition and light yield of a track in a paddle
 * Suitable geometry files require proper naming of the active volume (see CheckIfSensitive) and copy numbers.
 */

class R3BNeuland : public R3BDetector
{
  public:
    /** Default constructor */
    R3BNeuland();

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param trans   position
     *@param rot     rotation */
    R3BNeuland(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot = TGeoRotation());

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param combi   position + rotation */
    explicit R3BNeuland(const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Standard constructor.
     *@param nDP     number of double planes
     *@param trans   position
     *@param rot     rotation */
    R3BNeuland(int nDP, const TGeoTranslation& trans, const TGeoRotation& rot = TGeoRotation());

    /** Standard constructor.
     *@param nDP     number of double planes
     *@param combi   position + rotation */
    explicit R3BNeuland(int nDP, const TGeoCombiTrans& combi = TGeoCombiTrans());

    auto ProcessHits(FairVolume* /*v*/ = nullptr) -> bool override;

    void EndOfEvent() override;

    void Print(Option_t* /*unused*/ = "") const override;

    void Reset() override;

    auto CheckIfSensitive(std::string name) -> bool override;

    [[nodiscard]] auto GetCollection(Int_t /*iColl*/) const -> TClonesArray* override { return nullptr; }

    void Register() override;

  private:
    R3B::OutputVectorConnector<R3BNeulandPoint> neuland_points_{ "NeulandPoints" }; //!
    R3BNeulandGeoPar* neuland_geo_par_ = nullptr;                                   //!
    std::map<int, int> trackid_pid_map_;

    /** Track information to be stored until the track leaves the active volume. */
    int track_id_ = 0;
    int paddle_id_ = 0;
    TLorentzVector pos_in_{};
    TLorentzVector pos_out_{};
    TLorentzVector mom_in_{};
    TLorentzVector mom_out_{};
    double time_ = 0.;
    double length_ = 0.;
    double energy_loss_ = 0.;
    double light_yield_ = 0.;
    bool is_last_hit_done_ = false;
    int particle_id_ = 0;
    int partent_particle_id_ = 0;

    void ResetValues();

    void WriteParameterFile();

    ClassDefOverride(R3BNeuland, 3);
};

#endif // R3BNEULAND_H
