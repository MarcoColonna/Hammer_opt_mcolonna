#########################
#### Histograms demo ####
#########################

import sys
import time

from hammer import Hammer, IOBuffer, RecordType, BinContents

from math import sqrt

try:
    from ROOT.std import ifstream as ifstream
    from ROOT.std import ios as ios
except:
    try:
        from cppyy.gbl.std import ifstream as ifstream
        from cppyy.gbl.std import ios as ios
    except:
        pass

from cppyy.gbl.std import complex as cmplx
from cppyy.gbl.std import map, string

cd = cmplx['double']

def printHistogram(histo, keepErrors):
    for idx2 in range(0,5):
        for idx1 in range(0,6):
            print(f"{histo[idx1 * 5 + idx2].sumWi:10.4f}", end = "")
        print()    
        for idx1 in range(0,6):
            print(f"{histo[idx1 * 5 + idx2].n:10}", end = "")
        print()   
        if keepErrors:
            for idx1 in range(0,6):
                print(f"{sqrt(histo[idx1 * 5 + idx2].sumWi2):10.4f}", end = "")
            print("\n")    
    print()
    
if __name__ == '__main__':

    ##################################################
    #### Now reread histograms and play with them ####
    ##################################################

    ham = Hammer()
    buf = IOBuffer()
    
    inFile = ifstream("./DemoHistos.dat", ios.binary)
    ## Load from buffers (vs demo04.cc this is done just once, as a slightly different example)
    ham.setUnits("GeV")
    buf.load(inFile) ## First record saved was the run header
    if not ham.loadRunHeader(buf): 
        if (buf.start != nullptr):
            cppyy.ll.array_delete(buf.start)
        inFile.close()
        sys.exit(1)
    ham.initRun()
    buf.load(inFile) ## Next set of records saved were histograms
    while buf.kind == chr(RecordType.HISTOGRAM) or buf.kind == chr(RecordType.HISTOGRAM_DEFINITION):
        if buf.kind == chr(RecordType.HISTOGRAM):
            ham.loadHistogram(buf)
        else:
            ham.loadHistogramDefinition(buf)
        buf.load(inFile)
    inFile.close()
    ## Begin run          
    for idW in range(0,6):    
        rwgtstart = time.process_time() 
        val = idW * 0.2
        ham.setWilsonCoefficients("BtoCTauNu", ({"S_qLlL": cd(0,val), "T_qLlL": cd(val/4.,0)}))
        histoT = ham.getHistogram("Total Sum of Weights", "Scheme1");
        histoDs = ham.getHistogram("pEllVsQ2:D*", "Scheme1");
        histoD = ham.getHistogram("pEllVsQ2:D", "Scheme1")
        rwgtend = time.process_time()
        
        #histoT[0].sumWi
        print(f"\nReweighted histo to S_qLlL = {val:.1f}i, T_qLlL = {val/4.:.2f}: ")
        print("BD*TauNu: pEllVsQ2")
        printHistogram(histoDs, True)
        print("BDTauNu: pEllVsQ2")
        printHistogram(histoD, True)
        print(f"Total: {histoT[0].sumWi:.2f}")
        secsrwgt = rwgtend - rwgtstart
        print(f"Reweight Time: {secsrwgt:.5f}")