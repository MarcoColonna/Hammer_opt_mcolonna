///////////////////////
//// Weights demo ////
///////////////////////
#include <vector>
#include <chrono>
#include <fstream>

#include "Hammer/Hammer.hh"

using namespace std;

int main() { // NOLINT(bugprone-exception-escape)

    /////////////////////////////////////////////////////
    //// Now reread event weights and play with them ////
    /////////////////////////////////////////////////////

    Hammer::Hammer ham{};
    // If one wants to reprocess into a histogram
    // ham.addHistogram("Q2:D", {11}, false, {{0.,10.}});
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    for (size_t idW = 0; idW <= 5; ++idW) {
        auto rwgtstart = std::chrono::system_clock::now();
        ifstream inFile("./DemoWeights.dat", ios::binary);
        // Begin run
        double val = (static_cast<double>(idW)) * 0.2;
        ham.setUnits("GeV");
        ham.saveOptionCard("Opts02.yml", false);
        inFile >> buf;
        if (!ham.loadRunHeader(buf)) {
            delete[] buf.start;
            inFile.close();
            return EXIT_FAILURE;
        }
        ham.initRun();
        // Set the WCs
        ham.setWilsonCoefficients("BtoCTauNu", {{"S_qLlL", val*1i}, {"T_qLlL", val/4.}});
        // Note this is an incremental setting; only changes what you specify! One could do
        // ham.resetWilsonCoefficients("BtoCTauNu"); to reset to SM
        // ham.setWilsonCoefficients("BtoCTauNu", {1., val*1i, 0., 0., 0., val/4., 0., 0., 0., 0., 0.});
        // One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
        // Container for evtwgts
        vector<double> evtwgts;
        evtwgts.reserve(10000);
        inFile >> buf;
        size_t i = 0;
        while (buf.kind == Hammer::RecordType::EVENT) {
            if ((i + 1) % 1000 == 0) {
                cout << "." << flush;
            }
            ham.initEvent();
            ham.loadEventWeights(buf);
            // double evtwgt = ham.getWeight("Scheme1", "JM");
            double evtwgt = ham.getWeight("Scheme1");
            // We could have also done instead, without using evtIds:
            // double evtwgt = 1.;
            // auto wgtmap = ham.getWeights("Scheme1"); //get all process weights in the event map<HashId, double>
            // for(auto& elem : wgtmap){
            //     evtwgt *= elem.second;
            // }
            // Or we could have requested the weight restricted to specific subset of processes in the event
            // by passing the list of their IDs E.g. for two processes (proc1, proc2) one would write:
            // double evtwgt = ham.getWeight("Scheme1", {proc1, proc2});
            // Store the computed weight in a vector<double>, or do whatever you want with it!
            evtwgts.push_back(evtwgt);
            // For example, populate a histogram. NB the processEvent argument, that turns off weight recalculation.
            // ham.fillEventHistogram("Q2:D", {1.1});
            // ham.processEvent(Hammer::PAction::HISTOGRAMS);
            inFile >> buf;
            ++i;
        }
        cout << '\n'
             << "Reweighted " << evtwgts.size() << " events to S_qLlL = " << val << "i, T_qLlL = " << val / 4. << ": ";
        for (i = 0; i < 4; ++i) {
            cout << evtwgts[i] << ", ";
        }
        cout << evtwgts[4] << " ...." << '\n';
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        using float_seconds = std::chrono::duration<float>;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Reweight Time: " << secsrwgt.count() << '\n';

        //Next buffer element are rates
        while (buf.kind == Hammer::RecordType::RATE) {
            ham.loadRates(buf);
            inFile >> buf;
        }

        auto rateD = ham.getRate("B+D0barTau+Nu", "Scheme1");
        auto rateDs = ham.getRate("B+D*0barTau+Nu", "Scheme1");

        cout << "Reweighted B -> D Tau Nu rate: " << rateD << '\n';
        cout << "Reweighted B -> D* Tau Nu rate: " << rateDs << '\n';
        
        inFile.close();

        // Read out the populated histogram
        // auto shape = ham.getHistogramShape("Q2:D");
        // for (size_t idx1 = 0; idx1 < shape.size(); ++idx1){
        //     cout << shape[idx1] << "\t\t";
        // }
        // cout << '\n';
        // auto histo = ham.getHistogram("Q2:D", "Scheme1");
        // for (size_t idx1 = 0; idx1 < 10; ++idx1){
        //     cout << histo[idx1].n << "\t\t";
        // }
        // cout << '\n';
    }

} // int main()
