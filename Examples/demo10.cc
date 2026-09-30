//////////////////////////////
//// Histogram merge demo ////
//////////////////////////////
#include <vector>
#include <fstream>

#include "Hammer/Hammer.hh"

using namespace std;

//////////////////////////
//// Buffer functions ////
//////////////////////////

inline void printHistogram(Hammer::IOHistogram& histo, bool keepErrors = false) {
    for (size_t idx2 = 0; idx2 < 5; ++idx2){
            for (size_t idx1 = 0; idx1 < 6; ++idx1){
                cout << histo[idx1*5 + idx2].sumWi << "\t\t";
            }
            cout << '\n';
            for (size_t idx1 = 0; idx1 < 6; ++idx1){
                cout << histo[idx1*5 + idx2].n << "\t\t";
            }
            cout << '\n';
            if(keepErrors){
                for (size_t idx1 = 0; idx1 < 6; ++idx1){
                    cout << sqrt(histo[idx1*5 + idx2].sumWi2) << "\t\t";
                }
                cout << '\n' << '\n';
            }

        }
    cout << '\n';
}

int main() { // NOLINT(bugprone-exception-escape)

    //////////////////////////////////////////////
    //// Now reread both histograms and merge ////
    //////////////////////////////////////////////

    Hammer::Hammer ham{};
    vector<string> histofiles{"DemoHistos_1", "DemoHistos_2"};
    ham.setUnits("GeV");
    //First load and merge the headers
    for (auto& file : histofiles){
        Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
        ifstream inFile("./" + file + ".dat", ios::binary);
        inFile >> buf;
        if(!ham.loadRunHeader(buf, true)) { //Bool for merge set to true
            delete[] buf.start;
            inFile.close();
            return EXIT_FAILURE;
        }
        inFile.close();
    }
    ham.initRun();
    //Now load and merge the histograms
    for (auto& file : histofiles){
        Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
        ifstream inFile("./" + file + ".dat", ios::binary);
        inFile >> buf;
        inFile >> buf; //skip the header buffer entry
        while(buf.kind == Hammer::RecordType::HISTOGRAM || buf.kind == Hammer::RecordType::HISTOGRAM_DEFINITION) {
            if(buf.kind == Hammer::RecordType::HISTOGRAM) {
                ham.loadHistogram(buf, true);
            }
            else {
                ham.loadHistogramDefinition(buf, true);
            }
            inFile >> buf;
        }
        inFile.close();
    }


    auto histoDsMerge = ham.getHistogram("pEllVsQ2:D*", "Scheme1");

    cout << "Merged Histogram:" << '\n';
    printHistogram(histoDsMerge, true);

    //Now load original histo
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    ifstream inFile("./DemoHistos.dat", ios::binary);

    inFile >> buf;
    if(!ham.loadRunHeader(buf)) { //Bool for merge set to false (overwrites existing histo)
        delete[] buf.start;
        inFile.close();
        return EXIT_FAILURE;
    }

    inFile >> buf;
    while(buf.kind == Hammer::RecordType::HISTOGRAM || buf.kind == Hammer::RecordType::HISTOGRAM_DEFINITION) {
            if(buf.kind == Hammer::RecordType::HISTOGRAM) {
                ham.loadHistogram(buf); //Bool for merge set to false (overwrites existing histo)
            }
            else {
                ham.loadHistogramDefinition(buf);
            }
            inFile >> buf;
    }

    
    inFile.close();

    auto histoDsFull = ham.getHistogram("pEllVsQ2:D*", "Scheme1");
    cout << "Original Histogram:" << '\n';
    printHistogram(histoDsFull, true);


} // int main()
