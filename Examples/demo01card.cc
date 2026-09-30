///////////////////////
//// Weights demo ////
///////////////////////
#include "HepMCParser.hh"
#include "Hammer/Hammer.hh"

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
    ofstream outFile("./DemoWeights-Card.dat", ios::binary);
    Hammer::Hammer ham{};
    // Read inputs from cards
    ham.readCards("./data/Processes.yml", "./data/Options.yml");
    ham.addTotalSumOfWeights(); // adds "Total Sum of Weights" histo with auto bin filling
    ham.initRun();
    ham.saveOptionCard("Opts01-Card.yml", false);
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
    // save the total rates (optional)
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


} // int main()
