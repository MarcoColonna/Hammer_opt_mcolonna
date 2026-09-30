//////////////////////////////////////////////////////////
//// Histograms demo with FF eigenvectors and aliases ////
//////////////////////////////////////////////////////////
#include <vector>
#include <chrono>
#include <fstream>

#include "Hammer/Hammer.hh"

using namespace std;

//////////////////////////
//// Buffer functions ////
//////////////////////////

inline void printHistogram(Hammer::IOHistogram& histo, bool keepErrors = true) {
    for (size_t idx2 = 0; idx2 < 8; ++idx2){
            for (size_t idx1 = 0; idx1 < 7; ++idx1){
                cout << histo[idx1*8 + idx2].sumWi << "\t\t";
            }
            cout << '\n';
            for (size_t idx1 = 0; idx1 < 7; ++idx1){
                cout << histo[idx1*8 + idx2].n << "\t\t";
            }
            cout << '\n';
            if(keepErrors){
                for (size_t idx1 = 0; idx1 < 7; ++idx1){
                    cout << sqrt(histo[idx1*8 + idx2].sumWi2) << "\t\t";
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
        ifstream inFile("./DemoHistosFFAlias.dat", ios::binary);
        // Reload the histograms
        double val = (static_cast<double>(idW))*0.01;
        ham.setUnits("GeV");
        inFile >> buf;
        if(!ham.loadRunHeader(buf)) {
            delete[] buf.start;
            inFile.close();
            return EXIT_FAILURE;
        }
        ham.initRun();
        inFile >> buf;
        while(buf.kind == Hammer::RecordType::HISTOGRAM || buf.kind == Hammer::RecordType::HISTOGRAM_DEFINITION) {
            if(buf.kind == Hammer::RecordType::HISTOGRAM) {
                ham.loadHistogram(buf);
            }
            else {
                ham.loadHistogramDefinition(buf);
            }
            inFile >> buf;
        }
        auto rwgtstart = std::chrono::system_clock::now();
        //Set the FFs
        ham.setFFEigenvectors("BtoD", "BGLVar", {{"delta_ap1",val}, {"dir_2",val/3.}, {"delta_ap2",0.}});
        //Note this is an incremental setting; only changes what you specify! One could do
        //ham.resetFFEigenvectors("BtoD", "BGLVar"); to reset all deltas to zero
        //Cound also have done
        //ham.setFFEigenvectors("BtoD", "BGLVar",{{0.,val,val/3.,0.,0.,0.,0.,0.});
        //The output is row major flattened vector, so you need to know the dims
        auto histopEllQ2 = ham.getHistogram("pEllVsQ2:D", "SchemeVar");
        auto histoQ2 = ham.getHistogram("Q2:D", "SchemeVar");
        cout << '\n' << "Reweighted histo to delta_a1 = " << val << ", dir_2 = " << val/3. << ": " << '\n';
        cout << "BDEllNu: pEllVsQ2" << '\n';
        printHistogram(histopEllQ2);
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        using float_seconds = std::chrono::duration<float>;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Reweight Time: " << secsrwgt.count() << '\n';
        inFile.close();
    }

} // int main()
