/////////////////////////
//// Histograms demo ////
/////////////////////////
#include <chrono>

#include "HepMCParser.hh"
#include "Hammer/Tools/HammerRoot.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Tools/Pdg.hh"

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
// #include "TInterpreter.h"


using namespace std;


int main() { // NOLINT(bugprone-exception-escape)

    /////////////////////////////////////////////////////
    //// Now let's reprocess and add some histograms ////
    /////////////////////////////////////////////////////
    auto begin = std::chrono::system_clock::now();
    HepMC3::ReaderAsciiHepMC2 adapter{"./data/BCLepNu.hepmc"};
    if (adapter.failed()) {
        return EXIT_FAILURE;
    }
    HepMC3::Setup::set_print_warnings(false);
    // Outfile for Histos. More later!
    TFile* outFileH{new TFile{"./DemoHistos.root", "RECREATE"}};
    // prepares the tree
    TTree* tree{new TTree{"Hammer", "Hammer output Tree"}};
    tree->SetDirectory(outFileH);
    // gInterpreter->GenerateDictionary("std::set<size_t>","set");
    Hammer::RootIOBuffer treeBuffer;
    treeBuffer.maxLength = 64ul * 1024ul * 1024ul;
    treeBuffer.length = 0;
    treeBuffer.start = new UChar_t[static_cast<size_t>(treeBuffer.maxLength)];
    tree->Branch("type", &treeBuffer.kind, "type/B");
    tree->Branch("length", &treeBuffer.length, "length/I");
    tree->Branch("record", treeBuffer.start, "record[length]/b");

    Hammer::Hammer ham{};
    Hammer::IOBuffer fbBuffer;
    // Declare included processes. This time require a Tau -> Ell Nu Nu decay
    vector<string> includeDs = {"BD*TauNu", "TauEllNuNu"};
    vector<string> includeD = {"BDTauNu", "TauEllNuNu"};
    ham.includeDecay(includeDs);
    ham.includeDecay(includeD);
    ham.addFFScheme("Scheme1", {{"BD", "BLPR"}, {"BD*", "BLPR"}});
    ham.setFFInputScheme({{"BD", "ISGW2"}, {"BD*", "ISGW2"}});
    // Let's add a 2D histogram for the tau daughter lepton momentum versus q^2, 10 and 2 bins respectively
    ham.addHistogram("pEllVsQ2:D*", {6, 5}, true, {{0., 2.5}, {3., 12.}});
    ham.addHistogram("pEllVsQ2:D", {6, 5}, true, {{0., 2.5}, {3., 12.}});
    ham.addTotalSumOfWeights(); // adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV");
    ham.initRun();
    // Saves the FF scheme definitions & other run information
    fbBuffer = ham.saveRunHeader();
    treeBuffer = fbBuffer; // necessary to make sure buffer starts always at same location (alternatives?)
    tree->Fill();
    // Container for event Id for D and D* separately
    set<set<size_t>> evtIdDs;
    set<set<size_t>> evtIdD;
    size_t count = 0;
    // Loop over events
    HepMC3::GenEvent ge;
    for (size_t i = 0; i < 10000; ++i) {
        adapter.read_event(ge);
        if (adapter.failed()) {
            break;
        }
        if ((i + 1) % 1000 == 0) {
            cout << "processing event " << i + 1 << '\n';
        }
        auto processes = parseGenEvent(ge, {521, -521, 511, -511});
        if (!processes.empty()) {
            ham.initEvent();
            set<size_t> evtId;
            bool first = true;
            bool Dstar = true;
            for (auto& proc : processes) {
                auto procId = ham.addProcess(proc);
                if (procId == 0) { // Make sure process isn't forbidden
                    continue;
                }
                evtId.insert(proc.getId());
                if (!first) { // Only allow one tau decay to be binned.
                    continue;
                }
                // Build the process q^2 and pEll observable and determine the bin indices
                // Collect D* or D vertex particles and bin q2. Extract {Parent Particle, {Daughter Particles}}
                // This is done here using Hammer methods but in real life it will be done by the user using
                // other frameworks
                auto BParticles = proc.getParticlesByVertex("BD*TauNu");
                Dstar = true;
                string meson = "D*";
                if (BParticles.second.empty()) {
                    BParticles = proc.getParticlesByVertex("BDTauNu");
                    Dstar = false;
                    meson = "D";
                }
                double q2 = -1.;
                double pEll = -1.;
                for (auto& elem : BParticles.second) {
                    if (abs(elem.pdgId()) == Hammer::PID::DSTARPLUS || abs(elem.pdgId()) == Hammer::PID::DPLUS ||
                        abs(elem.pdgId()) == Hammer::PID::DSTAR || abs(elem.pdgId()) == Hammer::PID::D0) {
                        q2 = (BParticles.first.p() - elem.p()) * (BParticles.first.p() - elem.p());
                        break;
                    }
                }
                // Now get the tau decay vertex particles and bin pEll
                auto TauParticles = proc.getParticlesByVertex("TauEllNuNu");
                for (auto& elem2 : TauParticles.second) {
                    if (abs(elem2.pdgId()) == Hammer::PID::MUON || abs(elem2.pdgId()) == Hammer::PID::ELECTRON) {
                        pEll = elem2.p().p();
                        break;
                    }
                }
                // If everything exists, tell Hammer to add a bin entry for this event
                if (q2 >= 0. && pEll >= 0.) {
                    ham.fillEventHistogram("pEllVsQ2:" + meson, {pEll, q2});
                    first = false;
                }
            }
            if (!evtId.empty()) { // Process events with at least one tau
                ham.processEvent();
                if (Dstar) {
                    evtIdDs.insert(evtId);
                } else {
                    evtIdD.insert(evtId);
                }
                ++count;
            }
        }
    }
    // save the total rates (optional)
    fbBuffer = ham.saveRates();
    treeBuffer = fbBuffer;
    tree->Fill();
    // save the histograms
    for (auto& elem : ham.saveHistogram("Total Sum of Weights")) {
        treeBuffer = elem;
        tree->Fill();
    }
    for (auto& elem : ham.saveHistogram("pEllVsQ2:D*")) {
        treeBuffer = elem;
        tree->Fill();
    }
    for (auto& elem : ham.saveHistogram("pEllVsQ2:D")) {
        treeBuffer = elem;
        tree->Fill();
    }

    // write and close file
    tree->Write();
    outFileH->Write();
    outFileH->Close();

    auto end = std::chrono::system_clock::now();
    auto dur = end - begin;
    using float_seconds = std::chrono::duration<float>;
    auto secs = std::chrono::duration_cast<float_seconds>(dur);
    cout << "Histo Time: " << secs.count() << '\n';
    cout << "Events binned: " << count << '\n';
    cout << "B -> D* Histograms: " << evtIdDs.size() << '\n';
    cout << "B -> D Histograms: " << evtIdD.size() << '\n';

    // perform cleanup
    delete outFileH;
    delete[] treeBuffer.start;

} // int main()
