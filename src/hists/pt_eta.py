import awkward as ak
import hist
import numpy as np


pt_eta_reg_reg = (
                hist.Hist.new
                .Regular(100, 0, 100, name="pt_long", label="$p_T$ (GeV)")
                .Regular(50, 0, 2.5, name="eta_long", label=r"$\eta$")
                .IntCat([-10, 1, 10, 100, 1000], name="gen", label="Gen Tag")
                .IntCat([-1,1,10,100], name="qual", label="Quality Tag")
                #.IntCat([0,1,2,3,4,5], name="num_jet", label="NumJet")
                .Double()
            )


pt_eta_var_var_ele = (
                hist.Hist.new
                .Variable([2,3,4,5,7,10,20,45,75,1000], name="pt", label="Electron $p_T$ (GeV)")
                .Variable([0,0.8,1.4442,1.556,2.5], name="eta", label=r"Electron $\eta$")
                .IntCat([-10, 1, 10, 100, 1000], name="gen", label="Electron Gen Tag")
                .IntCat([-1,1,10,100], name="qual", label="Electron Quality Tag")
                #.IntCat([0,1,2,3,4,5], name="num_jet", label="NumJet")
                .Double()
            )


muon_hist_var_var_mu = (
                hist.Hist.new
                .Variable([2, 3, 4, 5, 7, 10, 20, 45, 75, 1000], name="pt", label="Muon $p_T$ (GeV)")
                .Variable([0, 0.9, 1.2, 2.1, 2.4], name="eta", label=r"Muon $\eta$")
                .IntCat([-10, 1, 10, 100, 1000], name="gen", label="Muon Gen Tag")
                .IntCat([-1,1,10,100], name="qual", label="Muon Quality Tag")
                #.IntCat([0,1,2,3,4,5], name="num_jet", label="NumJet")
                .Double()
            )