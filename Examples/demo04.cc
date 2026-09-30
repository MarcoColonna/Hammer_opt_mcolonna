/////////////////////////
//// Histograms demo ////
/////////////////////////
#include <vector>
#include <chrono>
#include <fstream>

#include "Hammer/Hammer.hh"

using namespace std;

//////////////////////////
//// Buffer functions ////
//////////////////////////

inline void printHistogram(Hammer::IOHistogram& histo, bool keepErrors = false) {
    for (size_t idx2 = 0; idx2 < 5; ++idx2) {
        for (size_t idx1 = 0; idx1 < 6; ++idx1) {
            cout << histo[idx1 * 5 + idx2].sumWi << "\t\t";
        }
        cout << '\n';
        for (size_t idx1 = 0; idx1 < 6; ++idx1) {
            cout << histo[idx1 * 5 + idx2].n << "\t\t";
        }
        cout << '\n';
        if (keepErrors) {
            for (size_t idx1 = 0; idx1 < 6; ++idx1) {
                cout << sqrt(histo[idx1 * 5 + idx2].sumWi2) << "\t\t";
            }
            cout << '\n' << '\n';
        }
    }
    cout << '\n';
}

int main() { // NOLINT(bugprone-exception-escape)

    //////////////////////////////////////////////////
    //// Now reread histograms and play with them ////
    //////////////////////////////////////////////////

    Hammer::Hammer ham{};
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    for (size_t idW = 0; idW <= 5; ++idW) {
        ifstream inFile("./DemoHistos.dat", ios::binary);
        // Reload the histograms
        double val = (static_cast<double>(idW)) * 0.2;
        ham.setUnits("GeV");
        inFile >> buf;
        if (!ham.loadRunHeader(buf)) {
            delete[] buf.start;
            inFile.close();
            return EXIT_FAILURE;
        }
        ham.initRun();
        inFile >> buf;
        while (buf.kind == Hammer::RecordType::HISTOGRAM || buf.kind == Hammer::RecordType::HISTOGRAM_DEFINITION) {
            if (buf.kind == Hammer::RecordType::HISTOGRAM) {
                ham.loadHistogram(buf);
            } else {
                ham.loadHistogramDefinition(buf);
            }
            inFile >> buf;
        }
        auto rwgtstart = std::chrono::system_clock::now();
        // Set the WCs
        ham.setWilsonCoefficients("BtoCTauNu", {{"S_qLlL", val*1i}, {"T_qLlL", val/4.}});
        // Note this is an incremental setting; only changes what you specify! One could do
        // ham.resetWilsonCoefficients("BtoCTauNu"); to reset to SM
        // One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
        // The output is row major flattened vector, so you need to know the dims, which can be obtained by getHistogramShape
        auto histoT = ham.getHistogram("Total Sum of Weights", "Scheme1");
        auto histoDs = ham.getHistogram("pEllVsQ2:D*", "Scheme1");
        auto histoD = ham.getHistogram("pEllVsQ2:D", "Scheme1");
        cout << '\n' << "Reweighted histo to S_qLlL = " << val << "i, T_qLlL = " << val / 4. << ": " << '\n';
        cout << "BD*TauNu: pEllVsQ2" << '\n';
        printHistogram(histoDs, true);
        cout << "BDTauNu: pEllVsQ2" << '\n';
        printHistogram(histoD, true);
        cout << "Total: " << histoT[0].sumWi << '\n';
        cout << "Total: " << sqrt(histoT[0].sumWi2) << '\n';
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        using float_seconds = std::chrono::duration<float>;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Load and Reweight Time: " << secsrwgt.count() << '\n';
        inFile.close();
    }

} // int main()
