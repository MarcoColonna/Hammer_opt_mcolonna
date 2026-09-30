import sys
import time

from hammer import Hammer
from hammer import PID

try:
    from ROOT.std import ios as ios
    from ROOT.std import ofstream as ofstream
except:
    try:
        from cppyy.gbl.std import ios as ios
        from cppyy.gbl.std import ofstream as ofstream
    except:
        pass
    
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
    
#########################
#### Histograms demo ####
#########################

if __name__ == '__main__':
    ###########################################
    #### Reprocess and add some histograms ####
    ###########################################

    begin = time.process_time()
    ## little hack to avoid linking errors in case HepMC has not been compiled with C++11 option on
    adapter = ReaderAsciiHepMC2("./data/BCLepNu.hepmc")
    if adapter.failed():
        sys.exit(1)
    ## Outfile for histos
    outFile = ofstream("./DemoHistos.dat", ios.binary)
    ham = Hammer()
    
    ## Declare included processes. This time require a Tau -> Ell Nu Nu decay
    includeDs = ["BD*TauNu", "TauEllNuNu"]
    includeD = ["BDTauNu", "TauEllNuNu"]
    ham.includeDecay(includeDs)
    ham.includeDecay(includeD)
    ## Declare FF schemes and input scheme
    ham.addFFScheme("Scheme1", [["BD", "BLPR"], ["BD*", "BLPR"]])
    ham.setFFInputScheme([["BD", "ISGW2"], ["BD*", "ISGW2"]])
    
    ## Now declare histograms
    ham.addHistogram("pEllVsQ2:D*", [6, 5], False, [[0., 2.5], [3., 12.]])
    ham.addHistogram("pEllVsQ2:D", [6, 5], False, [[0., 2.5], [3., 12.]])
    ## Turn on errors
    ham.keepErrorsInHistogram("pEllVsQ2:D*", True)
    ham.keepErrorsInHistogram("pEllVsQ2:D", True)
    ## Histogram compression
    ham.collapseProcessesInHistogram("pEllVsQ2:D*")
    ham.collapseProcessesInHistogram("pEllVsQ2:D")
    ##
    ham.addTotalSumOfWeights() ## adds "Total Sum of Weights" histo with auto bin filling
    ham.setUnits("GeV")
    ## For examples of WC specializations, see demo13.
    ham.initRun();
    ham.saveRunHeader().save(outFile)
    
    ## Container (set of frozensets) for event Id for D and D* separately 
    evtIdDs = set(())
    evtIdD = set(())
    count = 0
    ## Loop over events
    ge = hm.GenEvent()
    i = 0
    while i < 10000:
        adapter.read_event(ge)
        if adapter.failed():
            break
        if (i + 1) % 1000 == 0: 
            print(f"processing event {i+1}")
        ham.initEvent()
        processes = parseGenEvent(ge, [521, -521, 511, -511])
        evtId = set(())
        first = True
        Dstar = True
        for proc in processes:
            procId = ham.addProcess(proc)
            if procId == 0:
                continue
            evtId.add(proc.getId())     
            if not first: ## Only all one tau to be binned
                continue 
            ## Collect D* or D vertex particles. Extract {Parent Particle, {Daughter Particles}}
            BParticles = proc.getParticlesByVertex("BD*TauNu")
            Dstar = True
            meson = "D*"
            if len(BParticles[1]) == 0:
                BParticles = proc.getParticlesByVertex("BDTauNu")
                Dstar = False
                meson = "D"
            ## Now bin q2 observable
            q2bool = False
            pEllbool = False
            for elem in BParticles[1]:
                if abs(elem.pdgId()) == PID.DSTARPLUS or abs(elem.pdgId()) == PID.DPLUS or abs(elem.pdgId()) == PID.DSTAR or abs(elem.pdgId()) == PID.D0:
                    q2 = (BParticles[0].p() - elem.p()) * (BParticles[0].p() - elem.p())
                    q2bool = True
                    break
            ## Now do Tau decay vertex and pEll
            TauParticles = proc.getParticlesByVertex("TauEllNuNu")
            for elem2 in TauParticles[1]:
                if abs(elem2.pdgId()) == PID.MUON or abs(elem2.pdgId()) == PID.ELECTRON:
                    pEll = elem2.p().p()
                    pEllbool = True
                    break
            if q2bool and pEllbool:
                ham.fillEventHistogram(f"pEllVsQ2:{meson}", [pEll, q2])
                first = False
        if len(evtId) > 0:
            ham.processEvent()
            if Dstar:
                evtIdDs.add(frozenset(evtId))
            else:
                evtIdD.add(frozenset(evtId))
            count += 1    
        i += 1
    ham.saveHistogram("Total Sum of Weights").save(outFile)
    ham.saveHistogram("pEllVsQ2:D*").save(outFile)
    ham.saveHistogram("pEllVsQ2:D").save(outFile)
    outFile.close()
    init = time.process_time()
    secsinit = init - begin
    
    print(f"Histo Time: {secsinit:.5f}")
    print(f"Events binned: {count}")
    print(f"B -> D* Histograms: {len(evtIdDs)}")
    print(f"B -> D Histograms: {len(evtIdD)}")

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    