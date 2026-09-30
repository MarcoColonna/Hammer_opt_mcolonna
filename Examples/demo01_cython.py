######################
#### Weights demo ####
######################

import sys
import time

from hammer import Hammer, Process, Particle, FourMomentum

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
    ########################################################
    #### Initialize the sample and store tensor weights ####
    ########################################################

    begin = time.process_time()
    adapter = _ReaderAsciiHepMC2("./data/BCLepNu.hepmc")
    if adapter.failed():
        sys.exit(1)

    ham = Hammer()
    ## Declare included processes
    ham.include_decay(["BD*TauNu", "D*DPi"])
    ## Can be more restrictive, e.g. ["BD*TauNu", "D*DPi", "TauEllNuNu"]
    ham.include_decay("BDTauNu")
    ## Declare at least one FF scheme
    ham.add_ff_scheme("Scheme1", {"BD": "BLPR", "BD*": "BLPR"})
    ## Declare the input FF scheme
    ham.set_ff_input_scheme({"BD": "ISGW2", "BD*": "ISGW2"})
    ham.add_total_sum_of_weights()
    ham.set_units("GeV")
    ham.init_run()
    ham.save_option_card("./Opts01.yml", False)
    ham.save_header_card("./Card01.yml")

    with open("./DemoWeightsPY.dat", 'wb') as fout:
        ## Saves the FF scheme definitions & other run information
        ham.save_run_header().save(fout)
        ## Loop over events
        ge = _GenEvent()
        saved_events = 0
        i = 0
        while i < 10000:
            adapter.read_event(ge)
            if adapter.failed():
                break
            if (i + 1) % 1000 == 0:
                print(f"processing event {i+1}")
            ## Look for decay processes starting with B mesons
            processes = parse_gen_event(ge, [521, -521, 511, -511])
            if len(processes) > 0:
                has_processes = False
                ham.init_event()
                for elem in processes:
                    proc_id = ham.add_process(elem)
                    if proc_id != 0:
                        has_processes = True
                if has_processes:
                    saved_events += 1
                    ham.process_event()
                    ham.save_event_weights().save(fout)
            i += 1
        ## Save total rates (optional)
        ham.save_rates().save(fout)
        ## Save sum of weights of processed events
        ham.save_histogram("Total Sum of Weights").save(fout)

    elapsed = time.process_time() - begin
    print(f"Init Time: {elapsed:.5f}")
    print(f"Events Stored: {saved_events}")
    ham.save_references("./refs_demo01.bib")
