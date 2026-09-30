import sys
import time

from hammer import Hammer, BinSizes, BinRanges, Process, Particle, FourMomentum
import hammer.pdg as pdg

#########################
#### Histograms demo ####
#########################

# --- HepMC backend (cython-compatible: uses snake_case Hammer API) ---

try:
    from pyHepMC3 import HepMC3 as _hm
    _ReaderAsciiHepMC2 = _hm.ReaderAsciiHepMC2
    _hm.Setup.set_print_warnings(False)

    def _to_particle(p):
        m = p.momentum()
        return Particle(FourMomentum(m.e(), m.px(), m.py(), m.pz()), p.pdg_id())

    def _add_decay(proc, p, parent_id):
        daughters = []
        for po in p.children():
            daughters.append(proc.add_particle(_to_particle(po)))
            if po.status() == 2 and po.end_vertex() and po.pdg_id() != 111:
                _add_decay(proc, po, daughters[-1])
        proc.add_vertex(parent_id, daughters)

    def parse_gen_event(evt, pdg_codes):
        results = []
        for pi in evt.particles():
            if pi.status() == 2 and pi.pdg_id() in pdg_codes:
                p = Process()
                idx = p.add_particle(_to_particle(pi))
                _add_decay(p, pi, idx)
                results.append(p)
        return results

    _GenEvent = _hm.GenEvent

except ImportError:
    try:
        import pyhepmc as _hm
        from pyhepmc.io import ReaderAsciiHepMC2 as _ReaderAsciiHepMC2

        def _to_particle(p):
            m = p.momentum
            return Particle(FourMomentum(m.e, m.px, m.py, m.pz), p.pid)

        def _add_decay(proc, p, parent_id):
            daughters = []
            for po in p.children:
                daughters.append(proc.add_particle(_to_particle(po)))
                if po.status == 2 and po.end_vertex and po.pid != 111:
                    _add_decay(proc, po, daughters[-1])
            proc.add_vertex(parent_id, daughters)

        def parse_gen_event(evt, pdg_codes):
            results = []
            for pi in evt.particles:
                if pi.status == 2 and pi.pid in pdg_codes:
                    p = Process()
                    idx = p.add_particle(_to_particle(pi))
                    _add_decay(p, pi, idx)
                    results.append(p)
            return results

        _GenEvent = _hm.GenEvent

    except ImportError:
        from hammer.hepmc import HepMCParser as _HepMCParser

        class _ReaderAsciiHepMC2:
            def __init__(self, filename):
                self._parser = _HepMCParser(filename)
                self._failed = False

            def read_event(self, ge):
                event = self._parser.read_event()
                if event is None:
                    self._failed = True
                else:
                    ge._event = event

            def failed(self):
                return self._failed

        class _GenEvent:
            def __init__(self):
                self._event = None

        def _add_decay_hepmc(proc, pcode, parent_idx, evt):
            daughters = []
            for child_id in evt.particles[pcode].end_vertex(evt).children:
                child = evt.particles[child_id]
                idx = proc.add_particle(Particle(
                    FourMomentum(child.p.e, child.p.px, child.p.py, child.p.pz), child.pdg
                ))
                daughters.append(idx)
                if child.status_code == 2 and child.end_vertex_code != 0 and child.pdg != 111:
                    _add_decay_hepmc(proc, child_id, idx, evt)
            proc.add_vertex(parent_idx, daughters)

        def parse_gen_event(ge, pdg_codes):
            if ge._event is None:
                return []
            evt = ge._event
            results = []
            for pcode, pi in evt.particles.items():
                if pi.status_code == 2 and pi.pdg in pdg_codes:
                    p = Process()
                    parent_idx = p.add_particle(Particle(
                        FourMomentum(pi.p.e, pi.p.px, pi.p.py, pi.p.pz), pi.pdg
                    ))
                    _add_decay_hepmc(p, pcode, parent_idx, evt)
                    results.append(p)
            return results


if __name__ == '__main__':
    ###########################################
    #### Reprocess and add some histograms ####
    ###########################################

    begin = time.process_time()
    adapter = _ReaderAsciiHepMC2("./data/BCLepNu.hepmc")
    if adapter.failed():
        sys.exit(1)

    ham = Hammer()
    ## Declare included processes. This time require a Tau -> Ell Nu Nu decay
    ham.include_decay(["BD*TauNu", "TauEllNuNu"])
    ham.include_decay(["BDTauNu", "TauEllNuNu"])
    ## Declare FF schemes and input scheme
    ham.add_ff_scheme("Scheme1", {"BD": "BLPR", "BD*": "BLPR"})
    ham.set_ff_input_scheme({"BD": "ISGW2", "BD*": "ISGW2"})
    ## Now declare histograms
    ham.add_histogram("pEllVsQ2:D*", BinSizes([6, 5]), False, BinRanges([(0., 2.5), (3., 12.)]))
    ham.add_histogram("pEllVsQ2:D",  BinSizes([6, 5]), False, BinRanges([(0., 2.5), (3., 12.)]))
    ## Turn on errors
    ham.keep_errors_in_histogram("pEllVsQ2:D*", True)
    ham.keep_errors_in_histogram("pEllVsQ2:D", True)
    ## Histogram compression
    ham.collapse_processes_in_histogram("pEllVsQ2:D*")
    ham.collapse_processes_in_histogram("pEllVsQ2:D")
    ham.add_total_sum_of_weights()
    ham.set_units("GeV")
    ham.init_run()

    with open("./DemoHistosPY.dat", 'wb') as fout:
        ham.save_run_header().save(fout)

        ## Container (set of frozensets) for event Id for D and D* separately
        evtIdDs = set()
        evtIdD = set()
        count = 0
        ge = _GenEvent()
        i = 0
        while i < 10000:
            adapter.read_event(ge)
            if adapter.failed():
                break
            if (i + 1) % 1000 == 0:
                print(f"processing event {i+1}")
            ham.init_event()
            processes = parse_gen_event(ge, [521, -521, 511, -511])
            evtId = set()
            first = True
            Dstar = True
            for proc in processes:
                procId = ham.add_process(proc)
                if procId == 0:
                    continue
                evtId.add(proc.get_id())
                if not first:
                    continue
                ## Collect D* or D vertex particles. Extract {Parent Particle, {Daughter Particles}}
                BParticles = proc.get_particles_by_vertex("BD*TauNu")
                Dstar = True
                meson = "D*"
                if len(BParticles[1]) == 0:
                    BParticles = proc.get_particles_by_vertex("BDTauNu")
                    Dstar = False
                    meson = "D"
                ## Now bin q2 observable
                q2bool = False
                pEllbool = False
                for elem in BParticles[1]:
                    if abs(elem.pdg_id) in (pdg.DSTARPLUS, pdg.DSTAR, pdg.DPLUS, pdg.D0):
                        q2 = (BParticles[0].p - elem.p).mass2()
                        q2bool = True
                        break
                ## Now do Tau decay vertex and pEll
                TauParticles = proc.get_particles_by_vertex("TauEllNuNu")
                for elem2 in TauParticles[1]:
                    if abs(elem2.pdg_id) in (pdg.MUON, pdg.ELECTRON):
                        pEll = elem2.p.p()
                        pEllbool = True
                        break
                if q2bool and pEllbool:
                    ham.fill_event_histogram(f"pEllVsQ2:{meson}", [pEll, q2])
                    first = False
            if len(evtId) > 0:
                ham.process_event()
                if Dstar:
                    evtIdDs.add(frozenset(evtId))
                else:
                    evtIdD.add(frozenset(evtId))
                count += 1
            i += 1

        ham.save_histogram("Total Sum of Weights").save(fout)
        ham.save_histogram("pEllVsQ2:D*").save(fout)
        ham.save_histogram("pEllVsQ2:D").save(fout)

    init = time.process_time()
    secsinit = init - begin
    print(f"Histo Time: {secsinit:.5f}")
    print(f"Events binned: {count}")
    print(f"B -> D* Histograms: {len(evtIdDs)}")
    print(f"B -> D Histograms: {len(evtIdD)}")
