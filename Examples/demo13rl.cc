////////////////////////////////////////////////////////////
//// Weights reload and specialize into histograms demo ////
////////////////////////////////////////////////////////////
#include <chrono>
#include <fstream>

#include "HepMCParser.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Tools/Pdg.hh"

using namespace std;

int main() { // NOLINT(bugprone-exception-escape)

    /////////////////////////////////////////////////////
    //// Now reread event weights and play with them ////
    /////////////////////////////////////////////////////

    Hammer::Hammer ham{};

    // Load both precomputed weights, and the original hepMC (for kinematics)
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    HepMC3::ReaderAsciiHepMC2 adapter{"./data/BCLepNu.hepmc"};
    if (adapter.failed()) {
        return EXIT_FAILURE;
    }
    HepMC3::Setup::set_print_warnings(false);
    ifstream inFile("./DemoWeights.dat", ios::binary);

    // Begin run
    inFile >> buf;
    if (!ham.loadRunHeader(buf)) {
        delete[] buf.start;
        inFile.close();
        return EXIT_FAILURE;
    }
    // Reprocess into a histograms
    ham.addHistogram("pEll:D", {5}, false, {{0, 2.5}});
    ham.addHistogram("pEll:D*", {5}, false, {{0, 2.5}});
    ham.setUnits("GeV");
    ham.setOptions("ProcessCalc: {Rates: false}");
    // create a specialization
    ham.createWCSpecialization("RH", "BtoCTauNu", {"R2", "A"});
    ham.setWCSpecializationOrigin("RH", "BtoCTauNu", {{"SM", 1.0}});
    ham.setWCSpecializationBasis("RH", "BtoCTauNu",
                                 {{{"T_qLlL", 0.25}, {"S_qLlL", 1i}}, {{"V_qRlL", 1}, {"V_qLlL", -1}}});

    ham.applyWCSpecializationInHistogram("pEll:D", "RH");
    ham.applyWCSpecializationInHistogram("pEll:D*", "RH");
    ham.initRun();

    inFile >> buf;
    HepMC3::GenEvent ge;
    size_t count = 0;
    auto rwgtstart = std::chrono::system_clock::now();
    for (size_t i = 0; i < 10000; ++i) {
        adapter.read_event(ge);
        if (adapter.failed()) {
            break;
        }
        if ((i + 1) % 1000 == 0) {
            cout << "reprocessing event " << i + 1 << '\n';
        }
        ham.initEvent();
        // Reload processes for kinematics (2 per event)
        auto processes = parseGenEvent(ge, {521, -521, 511, -511});
        bool first = true;
        bool reprocess = false;
        double pEll = 0.;
        string meson;
        for (auto& proc : processes) {
            auto procId = ham.addProcess(proc);
            if (procId == 0) { // Make sure process isn't forbidden
                ham.removeProcess(procId);
                continue;
            }
            if (!first) { // Only allow one tau decay to be binned.
                ham.removeProcess(procId);
                continue;
            }
            // Collect D* or D vertex particles and bin pEll. Extract {Parent Particle, {Daughter Particles}}
            auto BParticles = proc.getParticlesByVertex("BD*TauNu");
            meson = "D*";
            if (BParticles.second.empty()) {
                BParticles = proc.getParticlesByVertex("BDTauNu");
                meson = "D";
            }
            // Now get the tau momentum
            for (auto& elem2 : BParticles.second) {
                if (abs(elem2.pdgId()) == Hammer::PID::TAU) {
                    pEll = elem2.p().p();
                    first = false;
                    reprocess = true;
                    break;
                }
            }
            ham.removeProcess(procId);
        }
        if (reprocess) {
            // Populate a histogram. NB the processEvent argument, that turns off weight recalculation.
            ham.loadEventWeights(buf);
            ham.fillEventHistogram("pEll:" + meson, {pEll});
            ham.processEvent(Hammer::PAction::HISTOGRAMS);
            inFile >> buf;
            ++count;
        }
    }
    cout << "Specialized and binned " << count << " events." << '\n';
    ;
    auto rwgtend = std::chrono::system_clock::now();
    auto durrwgt = rwgtend - rwgtstart;
    using float_seconds = std::chrono::duration<float>;
    auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
    cout << "Reprocess Time: " << secsrwgt.count() << '\n' << '\n';
    inFile.close();

    for (size_t idx = 0; idx < 5; ++idx) {
        double val = static_cast<double>(idx) * 0.2;
        ham.setWilsonCoefficients("BtoCTauNu@RH", {{"R2", val}, {"A", 0.5}});
        // Read out the populated histogram
        auto histoD = ham.getHistogram("pEll:D", "Scheme1", "RH");
        for (size_t idx1 = 0; idx1 < 5; ++idx1) {
            cout << histoD[idx1].sumWi << "\t\t";
        }
        cout << '\n';
        auto histoDs = ham.getHistogram("pEll:D*", "Scheme1", "RH");
        for (size_t idx1 = 0; idx1 < 5; ++idx1) {
            cout << histoDs[idx1].sumWi << "\t\t";
        }
        cout << '\n' << '\n';
    }
} // int main()
