import uproot
import matplotlib.pyplot as plt
import awkward as ak

# Open ROOT file and access tree
file = uproot.open("output.root")
tree = file["tout"]

# Load arrays
L1_DoubleMu0_Upt8_SQ_er2p0 = tree["L1_DoubleMu0_Upt8_SQ_er2p0"].array()
L1_DoubleMu0_Upt7_SQ_er2p0 = tree["L1_DoubleMu0_Upt7_SQ_er2p0"].array()
L1_DoubleMu_15_7 = tree["L1_DoubleMu_15_7"].array()
L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 = tree["L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7"].array()
L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18 = tree["L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18"].array()
L1_DoubleMu8_SQ = tree["L1_DoubleMu8_SQ"].array()
L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6 = tree["L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6"].array()
L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = tree["L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"].array()
L1_DoubleMu4p5_SQ_OS_dR_Max1p2 = tree["L1_DoubleMu4p5_SQ_OS_dR_Max1p2"].array()
L1_DoubleMu0_Upt15_Upt7 = tree["L1_DoubleMu0_Upt15_Upt7"].array()
L1_DoubleMu0_Upt6_IP_Min1_Upt4 = tree["L1_DoubleMu0_Upt6_IP_Min1_Upt4"].array()
L1_DoubleMu0_Upt6_SQ_er2p0 = tree["L1_DoubleMu0_Upt6_SQ_er2p0"].array()


PV_x = tree["PV_x"].array()
PV_y = tree["PV_y"].array()
PV_z = tree["PV_z"].array()
SV_xe = tree["SV_xe"].array()
SV_ye = tree["SV_ye"].array()
SV_ze = tree["SV_ze"].array()
SV_chi2 = tree["SV_chi2"].array()
SV_prob = tree["SV_prob"].array()
SV_chi2Ndof = tree["SV_chi2Ndof"].array()
SV_lxy = tree["SV_lxy"].array()
SV_l3d = tree["SV_l3d"].array()
Muon_vtxIdxs = tree["Muon_vtxIdxs"].array()
Muon_saHits = tree["Muon_saHits"].array()
Muon_saMatchedStats = tree["Muon_saMatchedStats"].array()
Muon_muHits = tree["Muon_muHits"].array()
Muon_muChambs = tree["Muon_muChambs"].array()
Muon_muCSCDT = tree["Muon_muCSCDT"].array()
Muon_muMatch = tree["Muon_muMatch"].array()
Muon_muMatchedStats = tree["Muon_muMatchedStats"].array()
Muon_muExpMatchedStats = tree["Muon_muExpMatchedStats"].array()
Muon_muMatchedRPC = tree["Muon_muMatchedRPC"].array()
Muon_pixHits = tree["Muon_pixHits"].array()
Muon_stripHits = tree["Muon_stripHits"].array()
Muon_pixLayers = tree["Muon_pixLayers"].array()
Muon_trkLayers = tree["Muon_trkLayers"].array()
Muon_pt = tree["Muon_pt"].array()
Muon_eta = tree["Muon_eta"].array()
Muon_phi = tree["Muon_phi"].array()
#Muon_chi2Ndof = tree["Muon_chi2Ndof"].array()

# Plots
# PV_x
plt.hist(PV_x, bins=50, range=(-0.1,0.1), histtype='step', label='PV_x')
plt.xlabel("PV_x")
plt.ylabel("Entries")
plt.title("PV_x")
plt.legend()
plt.grid(False)
plt.savefig("PV_x_plot.png")
plt.show()
plt.clf()

#Pt
Muon_pt = ak.flatten(tree["Muon_pt"].array())

plt.hist(Muon_pt, bins=50, range=(0,40), histtype='step', label='Muon_pt')
plt.xlabel("Pt_{mu} (GeV)")
plt.ylabel("Events")
plt.title("Pt_{mu} (GeV)")
plt.legend()
plt.grid(False)
plt.savefig("Pt_plot.png")
plt.show()
plt.clf()

#Lxy
SV_lxy = ak.flatten(tree["SV_lxy"].array())

plt.hist(SV_lxy, bins=50, range=(0,10), histtype='step', label='SV_lxy')
plt.xlabel("lxy (cm)")
plt.ylabel("Events")
plt.title("lxy (cm)")
plt.legend()
plt.grid(False)
plt.savefig("lxy_plot.png")
plt.show()
plt.clf()
