/////////////////////////
//// Histograms demo ////
/////////////////////////
#include <chrono>

#include "Hammer/Tools/HammerRoot.hh"
#include "Hammer/Hammer.hh"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

int main() { // NOLINT(bugprone-exception-escape)

    //////////////////////////////////////////////////
    //// Now reread histograms and play with them ////
    //////////////////////////////////////////////////

    Hammer::Hammer ham{};

    TFile* inFile{new TFile{"./DemoHistos.root", "READ"}};
    // prepares the tree
    if (inFile == nullptr) { return EXIT_FAILURE; }
    TTree* tree;
    inFile->GetObject("Hammer",tree);
    TBranch* brecord = nullptr;
    TBranch* btype = nullptr;
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, new uint8_t[16UL*1024*1024]};
    tree->SetBranchAddress("record",buf.start, &brecord);
    tree->SetBranchAddress("type",&buf.kind, &btype);
    Long64_t nrecords = tree->GetEntries();
    vector<unique_ptr<TH2D>> histosD;
    vector<unique_ptr<TH2D>> histosDs;
    for(size_t idW = 0; idW <= 5; ++idW) {
        string nameD = string("pEllVsQ2:D_") + to_string(idW);
        string nameDstar = string("pEllVsQ2:D*_") + to_string(idW);
        histosD.push_back(make_unique<TH2D>(nameD.c_str(), nameD.c_str(), 6, 0., 2.5, 5, 3., 12.));
        histosDs.push_back(make_unique<TH2D>(nameDstar.c_str(), nameDstar.c_str(), 6, 0., 2.5, 5, 3., 12.));
    }
    for (size_t idW = 0; idW <= 5; ++idW) {
        auto rwgtstart = std::chrono::system_clock::now();
        // Reload the histograms
        double val = (static_cast<double>(idW))*0.2;
        ham.setUnits("GeV");
        auto entry = tree->LoadTree(0);
        brecord->GetEntry(entry);
        btype->GetEntry(entry);
        if(!ham.loadRunHeader(buf)) {
            inFile->Close();
            delete inFile;
            return EXIT_FAILURE;
        }
        ham.initRun();
        for(Int_t i = 1 /* header and rates before histos */; i< nrecords; ++i) {
            entry = tree->LoadTree(i);
            brecord->GetEntry(entry);
            btype->GetEntry(entry);
            if(buf.kind == Hammer::RecordType::HISTOGRAM_DEFINITION) {
                ham.loadHistogramDefinition(buf);
            }
            else if(buf.kind == Hammer::RecordType::HISTOGRAM) {
                ham.loadHistogram(buf);
            }
        }
        //Set the WCs
        ham.setWilsonCoefficients("BtoCTauNu", {{"S_qLlL", val*1i}, {"T_qLlL", val/4.}});
        // One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
        ham.setHistogram2D("pEllVsQ2:D*", "Scheme1", *histosDs[idW]);
        ham.setHistogram2D("pEllVsQ2:D", "Scheme1", *histosD[idW]);
        auto histoT = ham.getHistogram("Total Sum of Weights", "Scheme1");
        cout << '\n' << "Reweighted histo to S_qLlL = " << val << "i, T_qLlL = " << val/4. << ": " << histoT[0].sumWi << '\n';
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        using float_seconds = std::chrono::duration<float>;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Reweight Time: " << secsrwgt.count() << '\n';
    }

    auto* cD = new TCanvas("cD","The D histograms",200,10,600,400);
    for(size_t i = 0; i< histosD.size(); ++i) {
        if(i==0) {
            histosD[i]->Draw("SURF1");
        } else {
            histosD[i]->Draw("SURF1 SAME");
            cD->Update();
        }
    }
    cD->Print("histoD.pdf", "pdf");
    auto* cDs = new TCanvas("cDs","The D* histograms",200,10,600,400);
    for(size_t i = 0; i< histosDs.size(); ++i) {
        if(i==0) {
            histosDs[i]->Draw("SURF1");
        } else {
            histosDs[i]->Draw("SURF1 SAME");
            cDs->Update();
        }
    }
    cDs->Print("histoDs.pdf", "pdf");

    histosDs.clear();
    histosD.clear();
    tree->ResetBranchAddresses();
    inFile->Close();

    delete inFile;
    delete[] buf.start;
    delete cD;
    delete cDs;

} // int main()
