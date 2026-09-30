///////////////////////
//// Weights demo ////
///////////////////////
#include <chrono>

#include "HepMCParser.hh"
#include "Hammer/Tools/HammerRoot.hh"
#include "Hammer/Hammer.hh"


#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
// #include "TInterpreter.h"

using namespace std;

int main() { // NOLINT(bugprone-exception-escape)

    ////////////////////////////////////////////////////////
    //// Initialize the sample and store tensor weights ////
    ////////////////////////////////////////////////////////

    auto begin = std::chrono::system_clock::now();
    HepMC3::ReaderAsciiHepMC2 adapter{"./data/BCLepNu.hepmc"};
    if (adapter.failed()) {
        return EXIT_FAILURE;
    }
    HepMC3::Setup::set_print_warnings(false);
    // Outfile
    TFile* outFile{new TFile{"./DemoWeights.root", "RECREATE"}};
    // prepares the tree
    TTree* tree{new TTree{"Hammer", "Hammer output Tree"}};
    tree->SetDirectory(outFile);
    // gInterpreter->GenerateDictionary("std::set<size_t>","set");
    Int_t eventNumber = 0;
    array<ULong_t, 2> processIds;
    Hammer::RootIOBuffer treeBuffer;
    treeBuffer.init(64UL * 1024 * 1024);
    tree->Branch("ev_number", &eventNumber, "ev_number/I");
    tree->Branch("proc_ids", processIds.data(), "proc_ids[2]/l");
    tree->Branch("type", &treeBuffer.kind, "type/B");
    tree->Branch("length", &treeBuffer.length, "length/I");
    tree->Branch("record", treeBuffer.start, "record[length]/b");

    Hammer::Hammer ham{};
    Hammer::IOBuffer fbBuffer;
    // Declare included processes
    vector<string> include = {"BD*TauNu", "D*DPi"};
    // Can be more restrictive, e.g. {"BD*TauNu","D*DPi", "TauEllNuNu"}
    ham.includeDecay(include);
    ham.includeDecay(string("BDTauNu"));
    // For moment neglect light leptons
    // ham.includeDecay(string("BD*EllNu"));
    // Declare forbidden processes here: Only relevant for processes already included in the 'includeDecay'. E.g.
    // ham.includeDecay(string("BD*TauNu"));
    // vector<string> forbid = {"BD*TauNu","D*DGamma"};
    // ham.forbidDecay(forbid);
    // Declare at least one FF scheme. Can mix different parametrizations for each process.
    ham.addFFScheme("Scheme1", {{"BD", "BLPR"}, {"BD*", "BLPR"}});
    // Can add more schemes eg: ham.addFFScheme("Scheme2", {{"BD", "CLN"}, {"BD*", "BGL"}});
    // Declare the input FF scheme
    ham.setFFInputScheme({{"BD", "ISGW2"}, {"BD*", "ISGW2"}});
    ham.addTotalSumOfWeights(true); // adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV");
    ham.initRun();
    // Saves the FF scheme definitions & other run information
    fbBuffer = ham.saveRunHeader();
    treeBuffer = fbBuffer; // necessary to make sure buffer starts always at same location (alternatives?)
    tree->Fill();
    // Loop over events. Can be more than one process per event if there are two taus, or zero!
    HepMC3::GenEvent ge;
    size_t savedEvents = 0ul;
    for (size_t i = 0; i < 10000; ++i) {
        adapter.read_event(ge);
        if (adapter.failed()) {
            break;
        }
        eventNumber = static_cast<Int_t>(i + 1);
        processIds.fill(0);
        if ((i + 1) % 1000 == 0) {
            cout << "processing event " << i + 1 << '\n';
        }
        // look for decay processes starting with B mesons
        auto processes = parseGenEvent(ge, {521, -521, 511, -511});
        // can also use Hammer PDG constants for increased readability
        // auto processes = parseGenEvent(ge, {PID::BPLUS, PID::BMINUS, PID::BZERO, -PID::BZERO});
        if (!processes.empty()) {
            bool hasProcesses = false;
            ham.initEvent();
            for (size_t idx = 0; idx < processes.size(); ++idx) {
                auto procId =
                    ham.addProcess(processes[idx]); // addProcess return an unique Id of that specific process
                                                    // that can be used to retrieve process-specific information later
                if (procId != 0) {                  // Make sure process isn't forbidden and add it to the procIds list
                    hasProcesses = true;
                    processIds[idx] = procId;
                }
            }
            if (hasProcesses) { // Process & store events with at least one tau
                ++savedEvents;
                ham.processEvent();
                fbBuffer = ham.saveEventWeights();
                treeBuffer = fbBuffer;
                tree->Fill();
            }
        }
    }
    // save the total rates (optional)
    processIds.fill(0);
    eventNumber = 0;
    fbBuffer = ham.saveRates();
    treeBuffer = fbBuffer;
    tree->Fill();
    // save the sum of weights of processed events (the ones for which processEvent() has been called)
    for (auto& elem : ham.saveHistogram("Total Sum of Weights")) {
        treeBuffer = elem;
        tree->Fill();
    }


    // write and close file
    tree->Write();
    outFile->Write();
    outFile->Close();

    auto init = std::chrono::system_clock::now();
    auto durinit = init - begin;
    using float_seconds = std::chrono::duration<float>;
    auto secsinit = std::chrono::duration_cast<float_seconds>(durinit);
    cout << "Init Time: " << secsinit.count() << '\n';
    cout << "Events Stored: " << savedEvents << '\n';

    // perform cleanup
    delete outFile;
    delete[] treeBuffer.start;

} // int main()
