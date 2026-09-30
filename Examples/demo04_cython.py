#########################
#### Histograms demo ####
#########################

import sys
import time

from math import sqrt
from hammer import Hammer, IOBuffer, RecordType


def print_histogram(histo, keep_errors):
    for idx2 in range(0, 5):
        for idx1 in range(0, 6):
            print(f"{histo[idx1 * 5 + idx2].sum_wi:10.4f}", end="")
        print()
        for idx1 in range(0, 6):
            print(f"{histo[idx1 * 5 + idx2].n:10}", end="")
        print()
        if keep_errors:
            for idx1 in range(0, 6):
                print(f"{sqrt(histo[idx1 * 5 + idx2].sum_wi2):10.4f}", end="")
            print("\n")
    print()


if __name__ == '__main__':

    ##################################################
    #### Now reread histograms and play with them ####
    ##################################################

    ham = Hammer()
    buf = IOBuffer()

    ham.set_units("GeV")
    with open("./DemoHistosPY.dat", 'rb', buffering=0) as fin:
        ## Load from buffers: first record is the run header
        if not buf.load(fin):
            sys.exit(1)
        if not ham.load_run_header(buf):
            sys.exit(1)
        ham.init_run()
        ## Next set of records are histograms
        if not buf.load(fin):
            sys.exit(1)
        while buf.kind == RecordType.HISTOGRAM or buf.kind == RecordType.HISTOGRAM_DEFINITION:
            if buf.kind == RecordType.HISTOGRAM:
                ham.load_histogram(buf)
            else:
                ham.load_histogram_definition(buf)
            if not buf.load(fin):
                break

    ## Reweight histograms for different Wilson coefficient values
    for idW in range(0, 6):
        rwgtstart = time.process_time()
        val = idW * 0.2
        ham.set_wilson_coefficients("BtoCTauNu", {"S_qLlL": val * 1j, "T_qLlL": val / 4.})
        histoT  = ham.get_histogram("Total Sum of Weights", "Scheme1")
        histoDs = ham.get_histogram("pEllVsQ2:D*", "Scheme1")
        histoD  = ham.get_histogram("pEllVsQ2:D",  "Scheme1")
        rwgtend = time.process_time()

        print(f"\nReweighted histo to S_qLlL = {val:.1f}i, T_qLlL = {val / 4.:.2f}: ")
        print("BD*TauNu: pEllVsQ2")
        print_histogram(histoDs, True)
        print("BDTauNu: pEllVsQ2")
        print_histogram(histoD, True)
        print(f"Total: {histoT[0].sum_wi:.2f}")
        secsrwgt = rwgtend - rwgtstart
        print(f"Reweight Time: {secsrwgt:.5f}")
