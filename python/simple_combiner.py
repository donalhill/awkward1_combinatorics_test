import uproot4 as uproot
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import awkward1 as ak

def calc_invariant_mass(p1, p2):
    masses = np.sqrt(2*p1["pt"]*p2["pt"]*(np.cosh(p1["eta"] - p2["eta"]) - np.cos(p1["phi"] - p2["phi"])))

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Roman']})
rc('text', usetex=True)

#Points to home dir of project
path = "./"

#Open ROOT file from Pythia Delphes production with FCCSW
file = uproot.open(f"{path}/output/rootfiles/FCCDelphesOutput_100events.root")
tree = file['events'] #Get TTree from file

#Get the charged tracks
tr = tree.arrays(filter_name="pfcharged.core*")
tr["p"] = np.sqrt(tr["pfcharged.core.p4.px"]**2 + tr["pfcharged.core.p4.py"]**2 + tr["pfcharged.core.p4.pz"]**2)
tr["pt"] = np.sqrt(tr["pfcharged.core.p4.px"]**2 + tr["pfcharged.core.p4.py"]**2)
tr['eta'] = np.log((tr['p'] + tr['pfcharged.core.p4.pz'])/(tr['p'] - tr['pfcharged.core.p4.pz']))/2
tr['phi'] = np.arctan2(tr['pfcharged.core.p4.py'], tr['pfcharged.core.p4.px'])

#Pions
pi_cut = abs(tr["pfcharged.core.pdgId"]) == 211
pi = tr[pi_cut]
#Number of pions in each event
pi_sum = ak.num(pi_cut)
#Keep events with 2 or more pions
pi = pi[pi_sum >= 2]

#Make pion pairs per event
pi_pairs = ak.combinations(pi, 2)
#pt of first pion pair from first event
#print(pi_pairs[0]["pt"][0])

pi1, pi2 = ak.unzip(pi_pairs)

m_pipi = calc_invariant_mass(pi1, pi2)

plt.hist(ak.flatten(m_pipi),bins=400)
plt.show()
