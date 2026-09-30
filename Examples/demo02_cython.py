########################
#### Reweight demo  ####
########################

import sys
import time

from hammer import Hammer, IOBuffer, RecordType


if __name__ == '__main__':
    #####################################################
    #### Reread event weights and reweight with WCs  ####
    #####################################################

    ham = Hammer()
    buf = IOBuffer()

    for idW in range(0, 5):
        rwgtstart = time.process_time()
        val = idW * 0.2

        with open("./DemoWeightsPY.dat", 'rb', buffering=0) as fin:
            ham.set_units("GeV")
            ham.save_option_card("./Opts02.yml", False)

            if not buf.load(fin):
                sys.exit(1)
            if not ham.load_run_header(buf):
                sys.exit(1)
            ham.init_run()

            ## Set Wilson coefficients
            ham.set_wilson_coefficients("BtoCTauNu", {"S_qLlL": val * 1j, "T_qLlL": val / 4.})

            evtwgts = []
            if not buf.load(fin):
                sys.exit(1)
            while buf.kind == RecordType.EVENT:
                if (len(evtwgts) + 1) % 1000 == 0:
                    print(".", end="", flush=True)
                ham.init_event()
                ham.load_event_weights(buf)
                evtwgts.append(ham.get_weight("Scheme1"))
                if not buf.load(fin):
                    break

        print(f"\nReweighted {len(evtwgts)} events to S_qLlL = {val:.1f}i, T_qLlL = {val/4.:.2f}: ", end="")
        for i in range(0, min(4, len(evtwgts))):
            print(f"{evtwgts[i]:.6f}, ", end="")
        if len(evtwgts) >= 5:
            print(f"{evtwgts[4]:.6f}, ....")
        else:
            print()

        rwgtend = time.process_time()
        print(f"Reweight Time: {rwgtend - rwgtstart:.5f}")
