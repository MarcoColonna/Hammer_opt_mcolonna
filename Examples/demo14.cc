////////////////////////////////////////////////////
//// Histograms demo demo with WC specialization////
////////////////////////////////////////////////////
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

    /////////////////////////////////////////////////////////////////////////////
    //// Now reread specialized histograms and compare with general versions ////
    /////////////////////////////////////////////////////////////////////////////

    Hammer::Hammer ham{};
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    vector<string> infiles{"DemoHistos", "DemoHistosSpec"};
    double val = 1;
    ham.setUnits("GeV");
    for(auto& infile : infiles){
        ifstream inFile("./"+infile+".dat", ios::binary);
        // Reload the headers and merge
        inFile >> buf;
        if (!ham.loadRunHeader(buf, true)) {
            delete[] buf.start;
            inFile.close();
            return EXIT_FAILURE;
        }
        inFile.close();
    }
    ham.initRun();
    // Reload the histograms and merge
    for(auto& infile : infiles){
        ifstream inFile("./"+infile+".dat", ios::binary);
        inFile >> buf;
        inFile >> buf; //skip the header buffer entry
        while (buf.kind == Hammer::RecordType::HISTOGRAM || buf.kind == Hammer::RecordType::HISTOGRAM_DEFINITION) {
            if (buf.kind == Hammer::RecordType::HISTOGRAM) {
                ham.loadHistogram(buf, true);
            } else {
                ham.loadHistogramDefinition(buf, true);
            }
            inFile >> buf;
        }
        inFile.close();
    }
    // Set the WCs
    ham.setWilsonCoefficients("BtoCTauNu@JM", {{"R2", val}});
    ham.setWilsonCoefficients("BtoCTauNu", {{"S_qLlL", val*1i}, {"T_qLlL", val/4.}});
    // Note this is an incremental setting; only changes what you specify! One could do
    // ham.resetWilsonCoefficients("BtoCTauNu"); to reset to SM
    // One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
    // The output is row major flattened vector, so you need to know the dims
    vector<string> specs{"", "JM"}; 
    using float_seconds = std::chrono::duration<float>;
    vector<float_seconds> times{};
    for(auto& spec: specs){
        cout << "Specialization: " << (spec.empty() ? "None" : spec);
        auto rwgtstart = std::chrono::system_clock::now();
        //auto histoT = ham.getHistogram("Total Sum of Weights", "Scheme1");
        auto histoDs = ham.getHistogram("pEllVsQ2:D*", "Scheme1", spec);
        auto histoD = ham.getHistogram("pEllVsQ2:D", "Scheme1", spec);
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        times.push_back(secsrwgt);
        if(spec.empty()){
            cout << '\n' << "Reweighted histo to S_qLlL = " << val << "i, T_qLlL = " << val / 4. << ": " << '\n';
        } else if (spec == "JM"){
            cout << '\n' << "Reweighted histo to R2 = " << val << ": " << '\n';
        }
        cout << "BD*TauNu: pEllVsQ2" << '\n';
        printHistogram(histoDs, true);
        cout << "BDTauNu: pEllVsQ2" << '\n';
        printHistogram(histoD, true);
        //cout << "Total: " << histoT[0].sumWi << '\n';
        //cout << "Total: " << sqrt(histoT[0].sumWi2) << '\n';
        cout << "Reweight Time: " << secsrwgt.count() << '\n';
        cout << '\n';
    }
    cout << "Specialization speed factor: " << times[0]/times[1] << '\n';
} // int main()
