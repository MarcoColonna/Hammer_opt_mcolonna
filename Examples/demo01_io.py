import sys
import time

from hammer import Hammer

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

######################
#### Weights demo ####
######################

if __name__ == '__main__':
    ########################################################
    #### Initialize the sample and store tensor weights ####
    ########################################################

    begin = time.process_time()
    ## little hack to avoid linking errors in case HepMC has not been compiled with C++11 option on
    adapter = ReaderAsciiHepMC2("./data/BCLepNu.hepmc")
    if adapter.failed():
        sys.exit(1)
    ham = Hammer()
    ## Declare included processes
    include = ["BD*TauNu", "D*DPi"]
    ## Can be more restrictive, e.g. ["BD*TauNu", "D*DPi", "TauEllNuNu"]
    ham.includeDecay(include)
    ham.includeDecay("BDTauNu")
    ## For moment neglect light leptons
    ## ham.includeDecay("BD*EllNu")
    ## Declare forbidden processes here: Only relevant for processes already included in the 'includeDecay'. E.g.
    ## ham.includeDecay("BD*TauNu")
    ## forbid = ["BD*TauNu","D*DGamma"]
    ## ham.forbidDecay(forbid)
    ## Declare at least one FF scheme. Can mix different parametrizations for each process.
    ham.addFFScheme("Scheme1", [["BD", "BLPR"], ["BD*", "BLPR"]])
    ## Can add more schemes eg: ham.addFFScheme("Scheme2", [["BD", "CLN"], ["BD*", "BGL"]])
    ## Can add duplications of params in different schemes, denoted by Parametrization_token
    ## eg:  ham.addFFScheme("Scheme1", [["BD", "CLN_1"], ["BD*", "BGL_1"]])
    ## and  ham.addFFScheme("Scheme2", [["BD", "CLN_2"], ["BD*", "BGL_2"]])
    ## Declare the input FF scheme
    ham.setFFInputScheme([["BD", "ISGW2"], ["BD*", "ISGW2"]])
    ham.addTotalSumOfWeights() ## adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV")
    ## Specialize WCs to specific model subspace. See demo13
    ##    ham.createWCSpecialization("JM", "BtoCTauNu", ["R2"])
    ##    ham.setWCSpecializationOrigin("JM", "BtoCTauNu", {{"SM", 1.0}})
    ##    ham.setWCSpecializationCoord("JM", "BtoCTauNu", "R2", {{"T_qLlL", 0.25}, {"S_qLlL", 1i}})
    ##    ham.applyWCSpecializationInWeights("JM")
    ##    ## ham.setOptions("Hammer: {CalcGeneralWeights:true}")
    ham.initRun()
    ## Turn off rate calculation
    ## ham.setOptions("ProcessCalc: {Rates: false}")
    ham.saveOptionCard("./Opts01.yml", False)
    ham.saveHeaderCard("./Card01.yml")

    ## Outfile
    with open("./DemoWeights.dat", 'wb') as outFile:
        ## Saves the FF scheme definitions & other run information
        ham.saveRunHeader().save(outFile)
        ## Loop over events. Can be more than one process per event if there are two taus, or zero!
        ge = hm.GenEvent()
        savedEvents = 0
        i = 0
        while i < 10000:
            adapter.read_event(ge)
            if adapter.failed():
                break
            if (i + 1) % 1000 == 0: 
                print(f"processing event {i+1}")
            ## look for decay processes starting with B mesons
            processes = parseGenEvent(ge, [521, -521, 511, -511])
            ## can also use Hammer PDG constants for increased readability
            ## processes = parseGenEvent(ge, [PID.BPLUS, PID.BMINUS, PID.BZERO, -PID.BZERO])
            if len(processes) > 0:
                hasProcesses = False
                ham.initEvent()
                for elem in processes:
                    procId = ham.addProcess(elem) ## addProcess return an unique Id of that specific process
                                                ## that can be used to retrieve process-specific information later
                    if procId != 0:            ## Make sure process isn't forbidden
                        hasProcesses = True
                if hasProcesses:  ## Process & store events with at least one tau
                    savedEvents +=1
                    ham.processEvent()
                    ham.saveEventWeights().save(outFile)
            i += 1
        ## save the total rates (optional)
        ham.saveRates().save(outFile)
        ## save the sum of weights of processed events (the ones for which processEvent() has been called)
        ham.saveHistogram("Total Sum of Weights").save(outFile)
        outFile.close()
        init = time.process_time()
        secsinit = init - begin
        print(f"Init Time: {secsinit}")
        print(f"Events Stored: {savedEvents}")

        ## Export used bib references
        ham.saveReferences("./refs_demo01.bib")

