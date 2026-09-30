///////////////////////
//// Weights demo ////
///////////////////////
#include <chrono>
#include <fstream>

#include "HepMCParser.hh"
#include "Hammer/Hammer.hh"

using namespace std;

int main() { // NOLINT(bugprone-exception-escape)

    ////////////////////////////////////////////////////////
    //// Initialize the sample and store tensor weights ////
    ////////////////////////////////////////////////////////

    auto begin = std::chrono::system_clock::now();
    // little hack to avoid linking errors in case HepMC has not been compiled with C++11 option on
    HepMC3::ReaderAsciiHepMC2 adapter{"./data/BCLepNu.hepmc"};
    if (adapter.failed()) {
        return EXIT_FAILURE;
    }
    HepMC3::Setup::set_print_warnings(false);
    // Outfile
    ofstream outFile("./DemoWeights.dat", ios::binary);
    Hammer::Hammer ham{};
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
    // Can add duplications of params in different schemes, denoted by Parametrization_token
    // eg:  ham.addFFScheme("Scheme1", {{"BD", "CLN_1"}, {"BD*", "BGL_1"}});
    // and  ham.addFFScheme("Scheme2", {{"BD", "CLN_2"}, {"BD*", "BGL_2"}});
    // Display a list of all known FF parametrizations for all known class prefixes, or for a specific one:
    // ham.showAvailableFFParams();
    // ham.showAvailableFFParams("BtoD*");
    // Declare the input FF scheme
    ham.setFFInputScheme({{"BD", "ISGW2"}, {"BD*", "ISGW2"}});
    ham.addTotalSumOfWeights(); // adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV");
    // Specialize WCs to specific model subspace. See demo13
    //    ham.createWCSpecialization("JM", "BtoCTauNu", {"R2"});
    //    ham.setWCSpecializationOrigin("JM", "BtoCTauNu", {{"SM", 1.0}});
    //    ham.setWCSpecializationCoord("JM", "BtoCTauNu", "R2", {{"T_qLlL", 0.25}, {"S_qLlL", 1i}});
    //    ham.applyWCSpecializationInWeights("JM");
    //    //ham.setOptions("Hammer: {CalcGeneralWeights: true}");
    //    //ham.setOptions("Hammer: {GeneralRates: true}");
    ham.initRun();
    // Turn off rate calculation
    // ham.setOptions("ProcessCalc: {Rates: false}");
    ham.saveOptionCard("./Opts01.yml", false);
    ham.saveHeaderCard("./Card01.yml");
    // Saves the FF scheme definitions & other run information
    outFile << ham.saveRunHeader();
    // Loop over events. Can be more than one process per event if there are two taus, or zero!
    HepMC3::GenEvent ge;
    size_t savedEvents = 0ul;
    for (size_t i = 0; i < 10000; ++i) {
        adapter.read_event(ge);
        if (adapter.failed()) {
            break;
        }
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
            for (auto& elem : processes) {
                auto procId = ham.addProcess(elem); // addProcess return an unique Id of that specific process
                                                    // that can be used to retrieve process-specific information later
                if (procId != 0) {                  // Make sure process isn't forbidden
                    hasProcesses = true;
                }
            }
            if (hasProcesses) { // Process & store events with at least one tau
                ++savedEvents;
                ham.processEvent();
                outFile << ham.saveEventWeights();
            }
        }
    }
    // save the total rates
    outFile << ham.saveRates();
    // save the sum of weights of processed events (the ones for which processEvent() has been called)
    outFile << ham.saveHistogram("Total Sum of Weights");
    outFile.close();
    auto init = std::chrono::system_clock::now();
    auto durinit = init - begin;
    using float_seconds = std::chrono::duration<float>;
    auto secsinit = std::chrono::duration_cast<float_seconds>(durinit);
    cout << "Init Time: " << secsinit.count() << '\n';
    cout << "Events Stored: " << savedEvents << '\n';

    // Export used bib references
    ham.saveReferences("./refs_demo01.bib");
} // int main()
