#include "Neutron2DCalibr.h"

#include "TBranch.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TColor.h"
#include "TCutG.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"

#include "R3BNeulandCluster.h"
#include "R3BNeulandNeutron2DPar.h"

namespace Neuland
{
    Neutron2DCalibr::Neutron2DCalibr() {}

    void Neutron2DCalibr::SetClusterFile(const UInt_t nNeutrons, const TString& filename)
    {
        fHists[nNeutrons] =
            new TH2D(TString::UItoa(nNeutrons, 10), TString::UItoa(nNeutrons, 10) + "n", 50, 0, 2000, 35, 0, 70);
        fHists.at(nNeutrons)->GetXaxis()->SetTitle("Total Energy [MeV]");
        fHists.at(nNeutrons)->GetYaxis()->SetTitle("Number of Clusters");

        TFile* file = new TFile(filename, "READ");
        TTree* tree = (TTree*)file->Get("cbmsim");
        TBranch* branch = tree->GetBranch("NeulandClusters");
        TClonesArray* clusters = new TClonesArray("R3BNeulandCluster");
        branch->SetAddress(&clusters);

        const Int_t nEntries = tree->GetEntries();
        for (Int_t ei = 0; ei < nEntries; ei++)
        {
            branch->GetEntry(ei);

            Double_t Etot = 0.;
            Int_t validClusters = 0;
            const Int_t nClusters = clusters->GetEntries();

            for (Int_t ci = 0; ci < nClusters; ci++)
            {
                R3BNeulandCluster* cluster = (R3BNeulandCluster*)clusters->At(ci);
                if (cluster->GetE() > 0.)
                {
                    Etot += cluster->GetE();
                    validClusters++;
                }
            }

            if (Etot > 0)
            {
                fHists[nNeutrons]->Fill(Etot, validClusters);
            }
            else
            {
                fHists[nNeutrons]->Fill(-1, -1);
            }
        }

        // Some delete action here
    }

    void Neutron2DCalibr::Optimize()
    {
        const UInt_t nVars = 3;

        ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
        min->SetMaxFunctionCalls(100000000);
        min->SetMaxIterations(10000000);
        min->SetTolerance(0.05);

        ROOT::Math::Functor f(
            [&](const Double_t* d)
            {
                return WastedEfficiency(d);
            },
            nVars);
        min->SetFunction(f);

        Double_t step[nVars] = { 0.001, 0.5, 0.5 };
        Double_t variable[nVars] = { 0.04, 10, 3 };
        Double_t lower[nVars] = { 0.001, 1, 3 };
        Double_t upper[nVars] = { 10, 100, 6 };
        min->SetLimitedVariable(0, "slope", variable[0], step[0], lower[0], upper[0]);
        min->SetLimitedVariable(1, "distance", variable[1], step[1], lower[1], upper[1]);
        min->SetLimitedVariable(2, "distance offset", variable[2], step[2], lower[2], upper[2]);

        min->Minimize();

        std::cout << "Neutron2DCalibr::Optimize done!" << std::endl;
    }

    TCutG* Neutron2DCalibr::GetCut(const UInt_t nNeutrons, const Double_t k, const Double_t k0, const Double_t m)
    {
        if (!fCuts[nNeutrons])
        {
            fCuts[nNeutrons] = new TCutG(TString::UItoa(nNeutrons, 10), 4);
            fCuts.at(nNeutrons)->SetVarX("Total Energy [MeV]");
            fCuts.at(nNeutrons)->SetVarY("Number of Clusters");
        }

        Double_t y0;
        if (nNeutrons != 0)
        {
            y0 = (nNeutrons - 1) * k + k0;
        }
        else
        {
            y0 = 0;
        }

        const Double_t y3 = nNeutrons * k + k0;
        const Double_t x1 = y0 / m;
        const Double_t x2 = y3 / m;
        // std::cout << "(0, " << y0 << ") (" << x1 << ", 0) (" << x2 << ", 0) ( 0, " << y3 << ")" << std::endl;

        TCutG* cut = fCuts[nNeutrons];
        cut->SetPoint(0, 0, y0);
        cut->SetPoint(1, x1, 0);
        cut->SetPoint(2, x2, 0);
        cut->SetPoint(3, 0, y3);
        return cut;
    }

    Double_t Neutron2DCalibr::WastedEfficiency(const Double_t* d)
    {
        Double_t m = d[0];
        Double_t k = d[1];
        Double_t k0 = d[2];
        GetCut(0, k, k0, m);

        Double_t wasted_efficiency = 0;
        for (auto& nh : fHists)
        {
            const UInt_t nNeutrons = nh.first;
            wasted_efficiency += 1. - ((Double_t)GetCut(nNeutrons, k, k0, m)->IntegralHist(nh.second) /
                                       (Double_t)nh.second->GetEntries());
        }
        return wasted_efficiency;
    }

    void Neutron2DCalibr::Draw(const TString& img) const
    {
        gStyle->SetOptStat(0);
        // reverse viridis
        Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000 };
        Double_t red[9] = { 246. / 255., 144. / 255., 74. / 255., 35. / 255., 28. / 255.,
                            33. / 255.,  43. / 255.,  51. / 255., 26. / 255. };
        Double_t green[9] = { 222. / 255., 200. / 255., 180. / 255., 150. / 255., 118. / 255.,
                              87. / 255.,  55. / 255.,  24. / 255.,  9. / 255. };
        Double_t blue[9] = { 0. / 255.,   35. / 255.,  72. / 255., 101. / 255., 112. / 255.,
                             114. / 255., 112. / 255., 96. / 255., 30. / 255. };
        TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);

        TCanvas* c = new TCanvas("Neutron2DCalibr", "Neuland Neutron2D Calibr", 1000, (fHists.size() + 1) / 2 * 500);
        c->Divide(2, (fHists.size() + 1) / 2);

        for (auto& nh : fHists)
        {
            c->cd(nh.first);
            nh.second->Draw("colz");
            fCuts.at(nh.first)->Draw("same");
        }

        if (!(img == ""))
        {
            c->Print(img);
        }

        c->Draw();
    }

    void Neutron2DCalibr::Print(std::ostream& out) const
    {
        out << "\t";
        for (const auto& nh : fHists)
        {
            out << nh.first << "n\t";
        }
        out << "Purity";
        out << std::endl;

        for (const auto& nc : fCuts)
        {
            const UInt_t nOut = nc.first;
            const TCutG* cut = nc.second;

            out << nOut << "n:\t";

            Double_t sum = 0.;
            for (const auto& nh : fHists)
            {
                sum += ((Double_t)cut->IntegralHist(nh.second) / (Double_t)nh.second->GetEntries());
                out << ((Double_t)cut->IntegralHist(nh.second) / (Double_t)nh.second->GetEntries()) << "\t";
            }
            if (fHists.find(nOut) != fHists.end())
            {
                out << (Double_t)cut->IntegralHist(fHists.at(nOut)) / (Double_t)(fHists.at(nOut)->GetEntries()) / sum;
            }
            out << std::endl;
        }
    }

    void Neutron2DCalibr::WriteParameterFile(const TString& parFile) const
    {
        FairRuntimeDb* rtdb = FairRuntimeDb::instance();

        FairParRootFileIo* io = new FairParRootFileIo(kTRUE);
        io->open(parFile);
        rtdb->setOutput(io);

        R3BNeulandNeutron2DPar* par = (R3BNeulandNeutron2DPar*)rtdb->getContainer("R3BNeulandNeutron2DPar");

        rtdb->addRun(1);
        par->SetNeutronCuts(fCuts);
        par->setChanged();
        rtdb->writeContainers();

        rtdb->saveOutput();
        rtdb->print();
    }
}; // namespace
