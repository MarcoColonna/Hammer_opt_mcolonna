///////////////////////
//// Phase space demo ////
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
    HepMC3::ReaderAsciiHepMC2 adapter{"./data/BCTauNu_PS.hepmc"};
    if (adapter.failed()) {
        return EXIT_FAILURE;
    }
    HepMC3::Setup::set_print_warnings(false);
    Hammer::Hammer ham{};
    // Declare included processes
    vector<string> include = {"BD*TauNu", "D*DPi"};
    ham.includeDecay(include);
    ham.includeDecay(string("BDTauNu"));
    include = {"BD*EllNu", "D*DPi"};
    ham.includeDecay(include);
    ham.includeDecay(string("BDEllNu"));
    ham.addFFScheme("Scheme1", {{"BD", "BLPR"}, {"BD*", "BLPR"}});
    ham.setFFInputScheme({{"BD", ""}, {"BD*", ""}}); //When a vertex will be declared as purePS, one can leave the param name empty.
    // Pure PS vertices in denominator for b to c, but *not* in Tau decays! Boolean arg is for numerator.
    // The D* decays are merged with parent, and must also be declared!
    set<string> purePS{"BDTauNu", "BD*TauNu", "D*DPi", "BDEllNu", "BD*EllNu"};
    //    ham.addPurePSVertices(purePS, Hammer::WTerm::NUMERATOR);
    ham.addPurePSVertices(purePS, Hammer::WTerm::DENOMINATOR);
    ham.addTotalSumOfWeights(); // adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV");
    ham.initRun();
    set<set<size_t>> evtIds;
    size_t count = 0;
    HepMC3::GenEvent ge;
    for (size_t i = 0; i < 1000; ++i) {
        adapter.read_event(ge);
        if (adapter.failed()) {
            break;
        }
        if ((i + 1) % 100 == 0) {
            cout << "processing event " << i + 1 << '\n';
        }
        ham.initEvent();
        auto processes = parseGenEvent(ge, {521, -521, 511, -511});
        set<size_t> evtId;
        for (auto& elem : processes) {
            auto procId = ham.addProcess(elem);
            if (procId != 0) {
                evtId.insert(procId);
            }
        }
        if (!evtId.empty()) {
            evtIds.insert(evtId);
            ++count;
            ham.processEvent();
        }
    }
    auto init = std::chrono::system_clock::now();
    auto durinit = init - begin;
    using float_seconds = std::chrono::duration<float>;
    auto secsinit = std::chrono::duration_cast<float_seconds>(durinit);
    cout << "Init Time: " << secsinit.count() << '\n';
    cout << "Events Binned: " << count << '\n';
    cout << "Histograms: " << evtIds.size() << '\n';

    for (size_t idW = 0; idW <= 5; ++idW) {
        auto rwgtstart = std::chrono::system_clock::now();
        double val = (static_cast<double>(idW)) * 0.2;
        ham.setWilsonCoefficients("BtoCTauNu", {{"S_qLlL", val}, {"T_qLlL", val / 4.}});
        ham.setWilsonCoefficients("BtoCMuNu", {{"S_qLlL", -val}, {"T_qLlL", -val / 4.}});
        auto histo = ham.getHistogram("Total Sum of Weights", "Scheme1");
        cout << '\n'
             << "Reweighted histo to S_qLlL = " << val << ", T_qL1L = " << val / 4. << ": " << histo[0].sumWi << '\n';
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Reweight Time: " << secsrwgt.count() << '\n';
    }
} // int main()
