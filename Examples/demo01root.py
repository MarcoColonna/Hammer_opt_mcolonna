######################/
#### Weights demo ####
######################/
import sys
import time
from array import array
from ctypes import addressof, c_uint8, c_uint32

import ROOT
from hammer import Hammer, RootIOBuffer

if True:
    import cppyy.ll

try:
    from pyHepMC3 import HepMC3 as hm
    ReaderAsciiHepMC2 = hm.ReaderAsciiHepMC2
    from HepMCParser_HepMC3 import parseGenEvent
    hm.Setup.set_print_warnings(False)
except ImportError:
    try:
        import pyhepmc as hm
        from pyhepmc.io import ReaderAsciiHepMC2 as ReaderAsciiHepMC2
        from HepMCParser_pyhepmc import parseGenEvent
    except ImportError:
        import HepMCParser_hepmc as hm
        ReaderAsciiHepMC2 = hm.ReaderAsciiHepMC2
        parseGenEvent = hm.parseGenEvent

if __name__ == '__main__':

    ########################################################
    #### Initialize the sample and store tensor weights ####
    ########################################################

    begin = time.process_time()
    adapter = ReaderAsciiHepMC2("./data/BCLepNu.hepmc")
    if adapter.failed():
        sys.exit(1)
    ## Outfile
    outFile = ROOT.TFile.Open("./DemoWeights.root", "RECREATE")
    ## prepares the tree
    tree = ROOT.TTree("Hammer", "Hammer output Tree")
    tree.SetDirectory(outFile)
    eventNumber = array('l',[0])
    processIds = array('L',[0,0])
    treeBuffer = RootIOBuffer()
    treeBuffer.init(64ul * 1024ul * 1024ul)
    tree.Branch("ev_number", eventNumber, "ev_number/I")
    tree.Branch("proc_ids", processIds, "proc_ids[2]/l")
    tree.Branch("type", ROOT.addressof(treeBuffer, 'kind'), "type/B")
    tree.Branch("length", ROOT.addressof(treeBuffer, 'length'), "length/i")
    tree.Branch("record", treeBuffer.start, "record[length]/b")

    ham = Hammer()
    ## Declare included processes
    include = ["BD*TauNu", "D*DPi"]
    ## Can be more restrictive, e.g. ["BD*TauNu","D*DPi", "TauEllNuNu"]
    ham.includeDecay(include)
    ham.includeDecay("BDTauNu")
    ## For moment neglect light leptons
    ## ham.includeDecay("BD*EllNu")
    ## Declare forbidden processes here: Only relevant for processes already included in the 'includeDecay'. E.g.
    ## ham.includeDecay("BD*TauNu")
    ## forbid = ["BD*TauNu","D*DGamma"]
    ## ham.forbidDecay(forbid);
    ## Declare at least one FF scheme. Can mix different parametrizations for each process.
    ham.addFFScheme("Scheme1", [["BD", "BLPR"], ["BD*", "BLPR"]])
    ## Can add more schemes eg: ham.addFFScheme("Scheme2", [["BD", "CLN"], ["BD*", "BGL"]])
    ## Declare the input FF scheme
    ham.setFFInputScheme([["BD", "ISGW2"], ["BD*", "ISGW2"]])
    ham.addTotalSumOfWeights(True) ## adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV")
    ham.initRun()
    ## Saves the FF scheme definitions & other run information
    fbBuffer = ham.saveRunHeader()
    treeBuffer.__assign__(fbBuffer) ## necessary to make sure buffer starts always at same location (alternatives?)
    tree.Fill()
    ## Loop over events. Can be more than one process per event if there are two taus, or zero!
    ge = hm.GenEvent()
    savedEvents = 0
    i = 0
    while i < 10000:
        adapter.read_event(ge)
        if adapter.failed():
            break
        eventNumber = i + 1
        processIds[0] = 0
        processIds[1] = 0
        if eventNumber % 1000 == 0: 
                print(f"processing event {eventNumber}")
        ## look for decay processes starting with B mesons
        processes = parseGenEvent(ge, {521, -521, 511, -511})
        ## can also use Hammer PDG constants for increased readability
        ## auto processes = parseGenEvent(ge, {PID::BPLUS, PID::BMINUS, PID::BZERO, -PID::BZERO});
        if len(processes) > 0:
            hasProcesses = False
            ham.initEvent()
            for idx, elem in enumerate(processes):
                procId = ham.addProcess(elem) ## addProcess return an unique Id of that specific process
                                                    ## that can be used to retrieve process-specific information later
                if procId != 0:            ## Make sure process isn't forbidden
                    hasProcesses = True
                    processIds[idx] = procId
            if hasProcesses: ## Process & store events with at least one tau
                savedEvents +=1
                ham.processEvent()
                fbBuffer = ham.saveEventWeights()
                treeBuffer.__assign__(fbBuffer)
                tree.Fill()
        i += 1
    ## save the total rates (optional)
    processIds[0] = 0
    processIds[1] = 0
    eventNumber = 0
    fbBuffer = ham.saveRates()
    treeBuffer.__assign__(fbBuffer)
    tree.Fill()
    ## save the sum of weights of processed events (the ones for which processEvent() has been called)
    for elem in ham.saveHistogram("Total Sum of Weights"):
        treeBuffer = elem
        tree.Fill()



    ## write and close file
    tree.Write()
    outFile.Write()
    outFile.Close()

    init = time.process_time()
    secsinit = init - begin
    print(f"Init Time: {secsinit}")
    print(f"Events Stored: {savedEvents}")

    ## perform cleanup
    # cppyy.ll.array_delete(treeBuffer.start)
