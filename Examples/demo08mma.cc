//////////////////////////////////////////////////
//// Histograms demo; invokable for Mathematica////
//// Takes arguments for FFs                   ////
//////////////////////////////////////////////////
#include <vector>
#include <fstream>
#include <iomanip>

#include "Hammer/Hammer.hh"

using namespace std;

//////////////////////////
//// Buffer functions ////
//////////////////////////

inline void printHistogram(Hammer::IOHistogram& histo) {
    for (size_t idx2 = 0; idx2 < 5; ++idx2){
            for (size_t idx1 = 0; idx1 < 6; ++idx1){
                cout << histo[idx1*2 + idx2].sumWi << "\t\t";
            }
            cout << '\n';
            for (size_t idx1 = 0; idx1 < 6; ++idx1){
                cout << histo[idx1*2 + idx2].n << "\t\t";
            }
            cout << '\n';
        }
    cout << '\n';
}


inline void tabHistogram2D(Hammer::IOHistogram& histo, vector<size_t>& dims, vector<pair<double,double>>& ranges, const string& filename){
    ofstream file;
    file.open("./data/"+filename);
    for(size_t id1 = 0; id1 < dims[0]; ++id1){
        double valx = ranges[0].first + (static_cast<double>(id1)+0.5)/static_cast<double>(dims[0])*(ranges[0].second - ranges[0].first);
        for(size_t id2 = 0; id2 < dims[1]; ++id2){
            double valy = ranges[1].first + (static_cast<double>(id2)+0.5)/static_cast<double>(dims[1])*(ranges[1].second - ranges[1].first);
            file << setprecision(5)<< valx << "\t" << valy << "\t" << histo[id1*dims[1] + id2].sumWi << '\n';
        }
    }
    file.close();
}

inline void tabHistogram1D(Hammer::IOHistogram& histo, size_t dim, pair<double,double>& ranges, const string& filename){
    ofstream file;
    file.open("./data/"+filename);
    for(size_t id1 = 0; id1 < dim; ++id1){
        double valx = ranges.first + (static_cast<double>(id1)+0.5)/static_cast<double>(dim)*(ranges.second - ranges.first);
        file << setprecision(5)<< valx << "\t" << histo[id1].sumWi << '\n';
    }
    file.close();
}

int main(int /*unused*/, char *argv[]) { // NOLINT(bugprone-exception-escape)

    //////////////////////////////////////////////////
    //// Now reread histograms and play with them ////
    //////////////////////////////////////////////////

    Hammer::Hammer ham{};
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    ifstream inFile("./DemoHistosFF.dat", ios::binary);
    // Reload the histograms
    double val = stod(argv[1]);
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
    //Set the WCs
    ham.setFFEigenvectors("BtoD", "BGLVar", {{"delta_ap1",val}, {"delta_ap2",val/3.}});
    //One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
    //The output is row major flattened vector, so you need to know the dims
    vector<size_t> dims{6,5};
    vector<pair<double,double>> ranges{make_pair(0.,2.5),make_pair(0.012,12.)};

    auto histopEllQ2 = ham.getHistogram("pEllVsQ2:D", "SchemeVar");
    auto histoQ2 = ham.getHistogram("Q2:D", "SchemeVar");

    tabHistogram2D(histopEllQ2, dims, ranges, "mmadatapEllQ2_D_FF.dat");
    tabHistogram1D(histoQ2, 10, ranges[1], "mmadataQ2_D_FF.dat");


    inFile.close();

} // int main()
