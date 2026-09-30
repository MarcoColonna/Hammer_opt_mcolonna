#######################
#### Weights demo ####
#######################

import sys
import time

from hammer import Hammer, IOBuffer, RecordType

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

if __name__ == '__main__':

    #####################################################
    #### Now reread event weights and play with them ####
    #####################################################

    ham = Hammer()
    ## If one wants to reprocess into a histogram
    ## ham.addHistogram("Q2:D", [11], False, [[0.,10.]]);
    buf = IOBuffer()
    for idW in range(0,5):
        rwgtstart = time.process_time()
        inFile = ifstream("./DemoWeights.dat", ios.binary)
        ## Begin run
        val = idW * 0.2
        ham.setUnits("GeV")
        ham.saveOptionCard("Opts02.yml", False)
        buf.load(inFile)
        if not ham.loadRunHeader(buf):
            if (buf.start != nullptr):
                cppyy.ll.array_delete(buf.start)
            inFile.close()
            sys.exit(1)
        ham.initRun()
        ## Set the WCs
        ham.setWilsonCoefficients("BtoCTauNu", ({"S_qLlL": cd(0,val), "T_qLlL": cd(val/4.,0)}))
        ## Note this is an incremental setting; only changes what you specify! One could do
        ## ham.resetWilsonCoefficients("BtoCTauNu") to reset to SM
        ## ham.setWilsonCoefficients("BtoCTauNu", [1., val*1j, 0., 0., 0., val/4., 0., 0., 0., 0., 0.])
        ## One could also do e.g. ham.setWilsonCoefficients("BtoCMuNu", ...); if one cared about muon mode NP, etc.
        ## Container for evtwgts
        evtwgts = []
        buf.load(inFile)
        i = 0
        while (buf.kind == chr(RecordType.EVENT)):
            if ((i + 1) % 1000 == 0):
                print(".", end="", flush=True)
            ham.initEvent()
            ham.loadEventWeights(buf)
            ## evtwgt = ham.getWeight("Scheme1", "JM")
            evtwgt = ham.getWeight("Scheme1")
            ## We could have also done instead, without using evtIds:
            ## evtwgt = 1.;
            ## wgtmap = ham.getWeights("Scheme1") ##get all process weights in the event map<HashId, double>
            ## for key, value in wgtmap:
            ##     evtwgt *= value
            ## Or we could have requested the weight restricted to specific subset of processes in the event
            ## by passing the list of their IDs E.g. for two processes (proc1, proc2) one would write:
            ## evtwgt = ham.getWeight("Scheme1", [proc1, proc2]);
            ## Store the computed weight in a list, or do whatever you want with it!
            evtwgts.append(evtwgt)
            ## For example, populate a histogram. NB the processEvent argument, that turns off weight recalculation.
            ## ham.fillEventHistogram("Q2:D", [1.1])
            ## ham.processEvent(PAction.HISTOGRAMS)
            buf.load(inFile)
            i +=1
        print(f"\nReweighted {len(evtwgts)} events to S_qLlL = {val:.1f}i, T_qLlL = {val/4.:.2f}: ", end="")
        for i in range(0,4):
            print(f"{evtwgts[i]:.6f}, ", end="")
        print(f"{evtwgts[4]:.6f}, ....")
        rwgtend = time.process_time()
        secsrwgt = rwgtend - rwgtstart
        print(f"Reweight Time: {secsrwgt:.5f}")
        inFile.close()

        ## Read out the populated histogram
        ## auto shape = ham.getHistogramShape("Q2:D");
        ## for (size_t idx1 = 0; idx1 < shape.size(); ++idx1){
        ##     cout << shape[idx1] << "\t\t";
        ## }
        ## cout << endl;
        ## auto histo = ham.getHistogram("Q2:D", "Scheme1");
        ## for (size_t idx1 = 0; idx1 < 10; ++idx1){
        ##     cout << histo[idx1].n << "\t\t";
        ## }
        ## cout << endl;

