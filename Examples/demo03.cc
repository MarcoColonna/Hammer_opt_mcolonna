/////////////////////////
//// Histograms demo ////
/////////////////////////
#include <chrono>
#include <fstream>

#include "HepMCParser.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Tools/Pdg.hh"

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
    ofstream outFileH("./DemoHistos.dat", ios::binary);
    Hammer::Hammer ham{};
    // Declare included processes. This time require a Tau -> Ell Nu Nu decay
    vector<string> includeDs = {"BD*TauNu", "TauEllNuNu"};
    vector<string> includeD = {"BDTauNu", "TauEllNuNu"};
    ham.includeDecay(includeDs);
    ham.includeDecay(includeD);
    ham.addFFScheme("Scheme1", {{"BD", "BLPR"}, {"BD*", "BLPR"}});
    ham.setFFInputScheme({{"BD", "ISGW2"}, {"BD*", "ISGW2"}});
    // Let's add a 2D histogram for the tau daughter lepton momentum versus q^2, 10 and 2 bins respectively
    //     ham.addHistogram("Q2:D*", {10},false);
    //     ham.addHistogram("Q2:D", {10},false);
    ham.addHistogram("pEllVsQ2:D*", {6, 5}, false, {{0., 2.5}, {3., 12.}});
    ham.addHistogram("pEllVsQ2:D", {6, 5}, false, {{0., 2.5}, {3., 12.}});
    // Turn on errors handling globally
    // ham.setOptions("Histos: {KeepErrors: true}");
    // Errors by histogram
    ham.keepErrorsInHistogram("pEllVsQ2:D*", true);
    ham.keepErrorsInHistogram("pEllVsQ2:D", true);
    // Histogram compression
    ham.collapseProcessesInHistogram("pEllVsQ2:D*");
    ham.collapseProcessesInHistogram("pEllVsQ2:D");
    ham.addTotalSumOfWeights(); // adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV");
    // For examples of WC specializations, see demo13.
    ham.initRun();
    outFileH << ham.saveRunHeader();
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
        ham.initEvent();
        auto processes = parseGenEvent(ge, {521, -521, 511, -511});
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
            // Collect D* or D vertex particles and bin q2. Extract {Parent Particle, {Daughter Particles}}
            auto BParticles = proc.getParticlesByVertex("BD*TauNu");
            Dstar = true;
            string meson = "D*";
            if (BParticles.second.empty()) {
                BParticles = proc.getParticlesByVertex("BDTauNu");
                Dstar = false;
                meson = "D";
            }
            double q2 = 0.;
            double pEll = 0.;
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
            // Can also use setEventHistogramBin("pEllVsQ2:"+meson, {pEllbin, q2bin})
            // if bins are directly computed
            if (q2 > 0. && pEll > 0.) {
                ham.fillEventHistogram("pEllVsQ2:" + meson, {pEll, q2});
                first = false;
                // q2bin = static_cast <uint16_t> (floor(5.*(q2-3.)/(12. - 3.)));
                // pEllbin = static_cast<uint16_t>(floor(6. * (pEll - 0.) / (2.5 - 0.)));
            }
        }
        if (!evtId.empty()) { // Process & store events with at least one tau
            ham.processEvent();
            if (Dstar) {
                evtIdDs.insert(evtId);
            } else {
                evtIdD.insert(evtId);
            }
            ++count;
        }
    }
    // outFileH << ham.saveRates();
    outFileH << ham.saveHistogram("Total Sum of Weights");
    //    outFileH << ham.saveHistogram("Q2:D*");
    //    outFileH << ham.saveHistogram("Q2:D");
    outFileH << ham.saveHistogram("pEllVsQ2:D*");
    outFileH << ham.saveHistogram("pEllVsQ2:D");
    outFileH.close();

    auto end = std::chrono::system_clock::now();
    auto dur = end - begin;
    using float_seconds = std::chrono::duration<float>;
    auto secs = std::chrono::duration_cast<float_seconds>(dur);
    cout << "Histo Time: " << secs.count() << '\n';
    cout << "Events binned: " << count << '\n';
    cout << "B -> D* Histograms: " << evtIdDs.size() << '\n';
    cout << "B -> D Histograms: " << evtIdD.size() << '\n';

} // int main()
