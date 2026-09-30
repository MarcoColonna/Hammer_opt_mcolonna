//////////////////////////////////////////////
//// Histograms demo with FF eigenvectors ////
//////////////////////////////////////////////
#include <chrono>
#include <fstream>

#include "HepMCParser.hh"
#include "Hammer/Hammer.hh"

using namespace std;


int main() { // NOLINT(bugprone-exception-escape)

    /////////////////////////////////////////
    //// Process and add some histograms ////
    /////////////////////////////////////////
    auto begin = std::chrono::system_clock::now();
    HepMC3::ReaderAsciiHepMC2 adapter{"./data/BDEllNu.hepmc"};
    if (adapter.failed()) {
        return EXIT_FAILURE;
    }
    HepMC3::Setup::set_print_warnings(false);
    // Outfile for Histos. More later!
    ofstream outFileH("./DemoHistosFF.dat", ios::binary);
    Hammer::Hammer ham{};
    ham.includeDecay(string("BDEllNu"));
    ham.addFFScheme("SchemeVar", {{"BD", "BGLVar"}});
    ham.setFFInputScheme({{"BD", "ISGW2"}});
    // Let's add a 2D histogram for lepton momentum versus q^2, 10 and 2 bins respectively
    ham.addHistogram("Q2:D", {12}, true, {{0., 12.}});
    ham.addHistogram("pEllVsQ2:D", {8, 7}, true, {{0., 2.5}, {0., 12.}});
    ham.addTotalSumOfWeights(); // adds "Total Sum of Weights" histo with auto bin filling
    // ham.keepErrorsInHistogram("pEllVsQ2:D", true);
    ham.setUnits("GeV");
    ham.initRun();
    outFileH << ham.saveRunHeader();
    // Container for event Id for D and D* separately
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
        for (auto& proc : processes) {
            auto procId = ham.addProcess(proc);
            if (procId == 0) { // Make sure process isn't forbidden
                continue;
            }
            evtId.insert(proc.getId());
            if (!first) { // Only allow one decay to be binned.
                continue;
            }
            // Collect D vertex particles and bin q2 and pEll. Extract {Parent Particle, {Daughter Particles}}
            auto BVer = proc.getParticlesByVertex("BDEllNu");

            auto pC = *find_if(BVer.second.begin(), BVer.second.end(),
                               [&](Hammer::Particle& p) -> bool { return (abs(p.pdgId()) >= 400); });
            double q2 = (BVer.first.p() - pC.p()) * (BVer.first.p() - pC.p());

            auto pL = *find_if(BVer.second.begin(), BVer.second.end(), [&](Hammer::Particle& p) -> bool {
                return ((abs(p.pdgId()) == 13) || (abs(p.pdgId()) == 11));
            });
            double pEll = pL.p().p();

            // If everythingin range, tell Hammer to add a bin entry for this event
            if (first) {
                ham.fillEventHistogram("pEllVsQ2:D", {pEll, q2});
                ham.fillEventHistogram("Q2:D", {q2});
                first = false;
            }
        }
        if (!evtId.empty()) { // Process & store events with at least one tau
            ham.processEvent();
            evtIdD.insert(evtId);
            ++count;
        }
    }
    // outFileH << ham.saveRates();
    outFileH << ham.saveHistogram("Total Sum of Weights");
    outFileH << ham.saveHistogram("Q2:D");
    outFileH << ham.saveHistogram("pEllVsQ2:D");
    outFileH.close();

    auto end = std::chrono::system_clock::now();
    auto dur = end - begin;
    using float_seconds = std::chrono::duration<float>;
    auto secs = std::chrono::duration_cast<float_seconds>(dur);
    cout << "Histo Time: " << secs.count() << '\n';
    cout << "Events binned: " << count << '\n';
    cout << "B -> D Histograms: " << evtIdD.size() << '\n';

} // int main()
