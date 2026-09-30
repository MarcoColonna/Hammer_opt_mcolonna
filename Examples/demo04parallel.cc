/////////////////////////
//// Histograms demo ////
/////////////////////////
#include <vector>
#include <chrono>
#include <fstream>
#include <thread>

#include "Hammer/Hammer.hh"


using namespace std;

static constexpr size_t NUM = 100;

inline void parallelFor(size_t numThreads, size_t start, size_t end, const std::function<void(size_t start, size_t end)>& work) {
    size_t workSize = (end - start) / numThreads;

    std::vector<std::thread> my_threads(numThreads);

    for (size_t i = 0; i < numThreads; ++i) {
        size_t batchStart = start + i * workSize;
        my_threads[i] = std::thread(work, batchStart, batchStart + workSize);
    }
    size_t lastStart = start + numThreads * workSize;
    work(lastStart, end);

    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
}

inline double calcHistoFun(Hammer::IOHistogram& histo) {
    double res = 0.;
    for (size_t idx2 = 0; idx2 < 5; ++idx2){
        for (size_t idx1 = 0; idx1 < 6; ++idx1){
            double err = sqrt(histo[idx1*5 + idx2].sumWi2);
            res += (err > 0.) ? sqrt(fabs(histo[idx1*5 + idx2].sumWi - static_cast<double>(histo[idx1*5 + idx2].n)))/err : 0.;
        }
    }
    return res;
}

using float_seconds = std::chrono::duration<float>;

int main() { // NOLINT(bugprone-exception-escape)

    //////////////////////////////////////////////////
    //// Now reread histograms and play with them ////
    //////////////////////////////////////////////////

    Hammer::Hammer ham{};
    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, nullptr};
    ifstream inFile("./DemoHistos.dat", ios::binary);
    // Reload the histograms
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
    inFile.close();
    array<double, NUM> serialResults;
    array<double, NUM> parallelResults;
    auto fun = [](size_t id, Hammer::Hammer* pHam, array<double, NUM>* results) -> void {
        double val = (static_cast<double>(id)) * 0.001 - 10.;
        // Set the WCs
        // "SM", "S_qLlL", "S_qRlL", "V_qLlL", "V_qRlL", "T_qLlL", "S_qLlR", "S_qRlR", "V_qLlR", "V_qRlR", "T_qRlR"
        // pHam->setWilsonCoefficientsLocal("BtoCTauNu", {1.,  val * 1i, 0., 0., 0., val / 4., 0., 0., 0., 0., 0.});
        pHam->setWilsonCoefficientsLocal("BtoCTauNu", {{"S_qLlL", val * 1i}, {"T_qLlL", val / 4.}});
        // Note this is an incremental setting; only changes what you specify! One could do
        // ham.resetWilsonCoefficients("BtoCTauNu"); to reset to SM
        // One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
        // The output is row major flattened vector, so you need to know the dims
        auto histoDs = pHam->getHistogram("pEllVsQ2:D*", "Scheme1");
        auto histoD = pHam->getHistogram("pEllVsQ2:D", "Scheme1");

        (*results)[id] = calcHistoFun(histoDs) + calcHistoFun(histoD);
    };
    cout << "Starting serial evaluation of " << NUM << " evaluations of two 2D histograms (2*30 bins w/errors, collapsed):" << '\n';
    auto rwgtstart = std::chrono::system_clock::now();
    for (size_t idW = 0; idW < NUM; ++idW) {
        fun(idW, &ham, &serialResults);
    }
    auto rwgtend = std::chrono::system_clock::now();
    auto durrwgt = rwgtend - rwgtstart;
    auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
    cout << "Serial Reweight Time: " << secsrwgt.count() << '\n';
    size_t tmp = std::thread::hardware_concurrency();
    for (size_t n = 2; n <= 4; n *= 2) {
        size_t numThreads = tmp == 0 ? n : min(n, tmp);
        cout << "Starting parallel evaluation of " << NUM
             << " evaluations of two 2D histograms (2*30 bins w/errors, collapsed) with " << numThreads
             << " threads:" << '\n';
        rwgtstart = std::chrono::system_clock::now();
        parallelFor(numThreads, 0, NUM, [&](size_t start, size_t end) {
            for (size_t idW = start; idW < end; ++idW) {
                fun(idW, &ham, &parallelResults);
            }
        });
        rwgtend = std::chrono::system_clock::now();
        durrwgt = rwgtend - rwgtstart;
        secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Parallel Reweight Time with " << numThreads << " threads: " << secsrwgt.count() << '\n';
    }
    cout << "Now comparing parallel vs serial results: ";
    double avgDiscrepancy = 0.;
    for (size_t i = 0; i < NUM; ++i) {
        double minval = min(serialResults[i], parallelResults[i]);
        double maxval = max(serialResults[i], parallelResults[i]);
        double term = fabs(serialResults[i]-parallelResults[i]);
        if(minval > 0.) {
            term /= (minval * NUM);
        }
        else if(maxval > 0.) {
            term /= (maxval * NUM);
        }
        else {
            term = 0.;
        }
        avgDiscrepancy += term;
    }
    cout << "average discrepancy = " << avgDiscrepancy << '\n';

} // int main()
