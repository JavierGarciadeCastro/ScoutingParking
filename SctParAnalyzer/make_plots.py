import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import mplhep
import re

#filename = "output_ctau-1-mA-2p00-mpi-10.root"
#filename = "output_ctau-10-mA-2p00-mpi-10.root"
filename = "output_ctau-100-mA-2p00-mpi-10.root"
#filename = "output_test.root"

file = uproot.open(filename)
tree = file["tout"]


# Extract parameters using regex
match = re.search(r"ctau-(\d+)-mA-(\d+p\d+)-mpi-(\d+)", filename)
if match:
    ctau = match.group(1)
    mA = match.group(2).replace("p", ".")
    mpi = match.group(3)
    title = f"ctau = {ctau} mm, mA = {mA} GeV, mpi = {mpi} GeV"
else:
    title = "Unknown Model Parameters"
    ctau = 'test'

print(title)


print('Reading the tree...')
###L1 seeds
#Single Muon
L1_SingleMu11_SQ14_BMTF = tree["L1_SingleMu11_SQ14_BMTF"].array()
L1_SingleMu10_SQ14_BMTF = tree["L1_SingleMu10_SQ14_BMTF"].array()

combined_L1_singleMu = [
  L1_SingleMu11_SQ14_BMTF,
  L1_SingleMu10_SQ14_BMTF
]

pass_L1_singleMu = [int(any(bits)) for bits in zip(*combined_L1_singleMu)]

#Double muon
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

combined_L1_doubleMu = [
    L1_DoubleMu0_Upt8_SQ_er2p0,
    L1_DoubleMu0_Upt7_SQ_er2p0,
    L1_DoubleMu_15_7,
    L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7,
    L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18,
    L1_DoubleMu8_SQ,
    L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6,
    L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4,
    L1_DoubleMu4p5_SQ_OS_dR_Max1p2,
    L1_DoubleMu0_Upt15_Upt7,
    L1_DoubleMu0_Upt6_IP_Min1_Upt4,
    L1_DoubleMu0_Upt6_SQ_er2p0
]
pass_L1_doubleMu = [int(any(bits)) for bits in zip(*combined_L1_doubleMu)]

combined_L1_all = [pass_L1_doubleMu, pass_L1_singleMu]


### HLT Scouting
#Single Muon
Scouting_SingleMu = tree["DST_PFScouting_SingleMuon_v*"].array()

combined_SingleMu_Scouting = [pass_L1_singleMu, Scouting_SingleMu]
pass_singleMu_Scouting = [int(all(bits)) for bits in zip(*combined_SingleMu_Scouting)]

#Double Muon
Scouting_DoubleMu = tree["DST_PFScouting_DoubleMuon_v*"].array()

combined_DoubleMu_Scouting = [pass_L1_doubleMu, Scouting_DoubleMu]
pass_DoubleMu_Scouting = [int(all(bits)) for bits in zip(*combined_DoubleMu_Scouting)]


### HLT Parking
#Single Muon
HLT_Mu9 = tree["HLT_Mu9_Barrel_L1HP10_IP6_v*"].array()
HLT_Mu10 = tree["HLT_Mu10_Barrel_L1HP11_IP6_v*"].array()
HLT_Mu7 = tree["HLT_Mu7_Barrel_L1HP8_IP6_v*"].array()
HLT_Mu8 = tree["HLT_Mu8_Barrel_L1HP9_IP6_v*"].array()
HLT_Mu0 = tree["HLT_Mu0_Barrel_L1HP6_IP6_v*"].array()
HLT_Mu6 = tree["HLT_Mu6_Barrel_L1HP7_IP6_v*"].array()

combined_PAT_HLT_singleMu = [HLT_Mu9, HLT_Mu10, HLT_Mu7, HLT_Mu8, HLT_Mu0, HLT_Mu6]

pass_PAT_HLT_singleMu = [int(any(bits)) for bits in zip(*combined_PAT_HLT_singleMu)]

combined_SingleMu_Parking = [pass_L1_singleMu, pass_PAT_HLT_singleMu]
pass_SingleMu_Parking = [int(all(bits)) for bits in zip(*combined_SingleMu_Parking)]

#Double Muon
HLT_DoubleMu4_3_LowMass_v12 = tree["HLT_DoubleMu4_3_LowMass_v*"].array()

combined_DoubleMu_Parking = [pass_L1_doubleMu, HLT_DoubleMu4_3_LowMass_v12]
pass_DoubleMu_Parking = [int(all(bits)) for bits in zip(*combined_DoubleMu_Parking)]

#Scouting variables
sct_PV_x = tree["sct_PV_x"].array()
sct_PV_y = tree["sct_PV_y"].array()
sct_PV_z = tree["sct_PV_z"].array()
sct_SV_xe = tree["sct_SV_xe"].array()
sct_SV_ye = tree["sct_SV_ye"].array()
sct_SV_ze = tree["sct_SV_ze"].array()
sct_SV_chi2 = tree["sct_SV_chi2"].array()
sct_SV_prob = tree["sct_SV_prob"].array()
sct_SV_chi2Ndof = tree["sct_SV_chi2Ndof"].array()
sct_SV_lxy = tree["sct_SV_lxy"].array()
sct_SV_l3d = tree["sct_SV_l3d"].array()
#sct_SV_selected = tree["sct_SV_selected"].array()

sct_Muon_pt = tree["sct_Muon_pt"].array()
sct_Muon_eta = tree["sct_Muon_eta"].array()
sct_Muon_phi = tree["sct_Muon_phi"].array()
sct_invMass = tree["sct_invMass"].array()
nsct_mu = tree["nsct_muons"].array()
sct_mu_selected = tree["sct_Muon_selected"].array()



#print(len(sct_SV_lxy))
#print(ak.to_list(sct_SV_lxy))
#print(ak.to_list(sct_mu_selected))
#print(selected_lxy)
#print(len(selected_lxy))
#print(sct_SV_selected)



##PAT variables:
PAT_Muon_pt = tree["PAT_Muon_pt"].array()
PAT_Muon_eta = tree["PAT_Muon_eta"].array()
PAT_Muon_phi = tree["PAT_Muon_phi"].array()
PAT_lxy = tree["PAT_lxy"].array()
PAT_Muon_vx = tree["PAT_vx"].array()
PAT_Muon_vy = tree["PAT_vy"].array()
PAT_invMass = tree["PAT_invMass"].array()
nPAT_mu = tree["nPAT_muons"].array()
PAT_mu_selected = tree["PAT_Muon_selected"].array()
PAT_pair_index = tree["PAT_pair_index"].array()
PAT_selected_index = tree["PAT_selected_index"].array()


##GEN variables:
Gen_pt = tree["Gen_pt"].array()
Gen_eta = tree["Gen_eta"].array()
Gen_phi = tree["Gen_phi"].array()
Gen_m = tree["Gen_m"].array()
Gen_vx = tree["Gen_vx"].array()
Gen_vy = tree["Gen_vy"].array()
Gen_vz = tree["Gen_vz"].array()
Gen_lxy = tree["Gen_lxy"].array()
Gen_status = tree["Gen_status"].array()
Gen_index = tree["Gen_index"].array()
Gen_pdgId = tree["Gen_pdgId"].array()
Gen_motherIndex = tree["Gen_motherIndex"].array()
Gen_motherPdgId = tree["Gen_motherPdgId"].array()
Gen_ct = tree["Gen_ct"].array()
ngenMuons = tree["Number_GenMuons"].array()

#print(ngenMuons)


Draw_plots = True
if Draw_plots:

  #################
  #Denominators
  #################
  
  print('Computing denominators ...')
  lxy_sct_with_2genMu = [x for x, nmu in zip(sct_SV_lxy, ngenMuons) if nmu >= 2]
  lxy_PAT_with_2genMu = [x for x, nmu in zip(PAT_lxy, ngenMuons) if nmu >= 2]
  
  pt_sct_with_2genMu = [x for x, nmu in zip(sct_Muon_pt, ngenMuons) if nmu >= 2]
  pt_PAT_with_2genMu = [x for x, nmu in zip(PAT_Muon_pt, ngenMuons) if nmu >= 2]
  
  print('Computing numerators ...')
  
  #################
  #Numerators
  #################
  
  #Events where there are 2 gen muons and we reconstruct them offline
  sct_selected_lxy = []

  for lxy_row, mu_sel in zip(sct_SV_lxy, sct_mu_selected):
    n_selected = sum(mu_sel)
    if n_selected >= 4 and len(lxy_row) >= 2:
      sct_selected_lxy.append(lxy_row[0])
      sct_selected_lxy.append(lxy_row[1])
    elif n_selected >= 2 and len(lxy_row) >= 1:
      sct_selected_lxy.append(lxy_row[0])

  lxy_sct_with_2genMu_and_2sct_mu = sct_selected_lxy
  lxy_PAT_with_2genMu_and_2PAT_mu = [[row[i] for i in index]for row, index in zip(PAT_lxy, PAT_pair_index)]

  pt_sct_with_2genMu_and_2sct_mu = [[pt for pt, keep in zip(row, mask) if keep]for row, mask in zip(sct_Muon_pt, sct_mu_selected)]
  pt_PAT_with_2genMu_and_2PAT_mu = [[row[i] for i in index] for row, index in zip(PAT_Muon_pt, PAT_selected_index)]
  
  print('Computing events that pass the triggers ...')

  lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers = []
  pt_sct_with_2genMu_and_2sct_mu_and_passTriggers = []
  
  lxy_PAT_with_2genMu_and_2PAT_mu_and_passTriggers = []
  pt_PAT_with_2genMu_and_2PAT_mu_and_passTriggers = []
  
  #Scouting
  for lxy_row, mu_sel, pt_row, flag in zip(sct_SV_lxy, sct_mu_selected, sct_Muon_pt, pass_DoubleMu_Scouting):
    n_selected = sum(mu_sel)
    if flag == 1:
      if n_selected >= 4 and len(lxy_row) >= 2:
        lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers.extend([lxy_row[0], lxy_row[1]])
      elif n_selected >= 2 and len(lxy_row) >= 1:
        lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers.append(lxy_row[0])
      pt_sct_with_2genMu_and_2sct_mu_and_passTriggers.append([pt for pt, keep in zip(pt_row, mu_sel) if keep])
  
  #PAT
  for lxy_row, idx_pair, pt_row, idx_sel, flag in zip(PAT_lxy, PAT_pair_index, PAT_Muon_pt, PAT_selected_index, pass_DoubleMu_Parking):
    if flag == 1:
        lxy_PAT_with_2genMu_and_2PAT_mu_and_passTriggers.append([lxy_row[i] for i in idx_pair])
        pt_PAT_with_2genMu_and_2PAT_mu_and_passTriggers.append([pt_row[i] for i in idx_sel])



  '''
  #Events where there are 2 gen muons, we reconstruct them offline and they pass the Double muon triggers
  lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers = [x for x, nmu, nsct, flag in zip(sct_SV_lxy, ngenMuons, nsct_mu, pass_DoubleMu_Scouting) if nmu >= 2 and nsct >=2 and flag ==1] 
  lxy_PAT_with_2genMu_and_2PAT_mu_and_passTriggers = [x for x, nmu, nPAT, flag in zip(PAT_lxy, ngenMuons, nPAT_mu, pass_DoubleMu_Parking) if nmu >= 2 and nPAT >=2 and flag ==1]
  
  pt_sct_with_2genMu_and_2sct_mu_and_passTriggers = [x for x, nmu, nsct, flag in zip(sct_Muon_pt, ngenMuons, nsct_mu, pass_DoubleMu_Scouting) if nmu >= 2 and nsct >=2 and flag ==1]
  pt_PAT_with_2genMu_and_2PAT_mu_and_passTriggers = [x for x, nmu, nPAT, flag in zip(PAT_Muon_pt, ngenMuons, nPAT_mu, pass_DoubleMu_Parking) if nmu >= 2 and nPAT >=2 and flag ==1]
  '''
  print('Flattening arrays ...')
  
  lxy_sct_with_2genMu = ak.flatten(lxy_sct_with_2genMu)
  lxy_PAT_with_2genMu = ak.flatten(lxy_PAT_with_2genMu)
  #lxy_sct_with_2genMu_and_2sct_mu = ak.flatten(lxy_sct_with_2genMu_and_2sct_mu) #Already flat
  lxy_PAT_with_2genMu_and_2PAT_mu = ak.flatten(lxy_PAT_with_2genMu_and_2PAT_mu)
  #lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers = ak.flatten(lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers) #Already flat
  lxy_PAT_with_2genMu_and_2PAT_mu_and_passTriggers = ak.flatten(lxy_PAT_with_2genMu_and_2PAT_mu_and_passTriggers)
  
  pt_sct_with_2genMu = ak.flatten(pt_sct_with_2genMu)
  pt_PAT_with_2genMu = ak.flatten(pt_PAT_with_2genMu)
  pt_sct_with_2genMu_and_2sct_mu = ak.flatten(pt_sct_with_2genMu_and_2sct_mu)
  pt_PAT_with_2genMu_and_2PAT_mu = ak.flatten(pt_PAT_with_2genMu_and_2PAT_mu)
  pt_sct_with_2genMu_and_2sct_mu_and_passTriggers = ak.flatten(pt_sct_with_2genMu_and_2sct_mu_and_passTriggers)
  pt_PAT_with_2genMu_and_2PAT_mu_and_passTriggers = ak.flatten(pt_PAT_with_2genMu_and_2PAT_mu_and_passTriggers)
  
  
  print('Plotting ...')

  nbins = 100
  plt.style.use(mplhep.style.CMS)

  ############################
  #Lxy Efficiency
  ############################

  bins_lxy = np.linspace(0,40,nbins + 1)
  # Denominator
  lxy_counts_sct_denom, _ = np.histogram(lxy_sct_with_2genMu, bins=bins_lxy)
  lxy_counts_PAT_denom, _ = np.histogram(lxy_PAT_with_2genMu, bins=bins_lxy)
  
  # Numerator
  lxy_counts_sct_num1, _ = np.histogram(lxy_sct_with_2genMu_and_2sct_mu, bins=bins_lxy)
  lxy_counts_PAT_num1, _ = np.histogram(lxy_PAT_with_2genMu_and_2PAT_mu, bins=bins_lxy)

  lxy_counts_sct_num2, _ = np.histogram(lxy_sct_with_2genMu_and_2sct_mu_and_passTriggers, bins=bins_lxy)
  lxy_counts_PAT_num2, _ = np.histogram(lxy_PAT_with_2genMu_and_2PAT_mu_and_passTriggers, bins=bins_lxy)


  lxy_efficiency1_sct = lxy_counts_sct_num1 / lxy_counts_sct_denom
  lxy_efficiency1_PAT = lxy_counts_PAT_num1 / lxy_counts_PAT_denom

  lxy_efficiency2_sct = lxy_counts_sct_num2 / lxy_counts_sct_denom
  lxy_efficiency2_PAT = lxy_counts_PAT_num2 / lxy_counts_PAT_denom


  bin_centers_lxy = 0.5 * (bins_lxy[:-1] + bins_lxy[1:])

  plt.step(bin_centers_lxy, lxy_efficiency1_sct, where='mid', label='Scouting RECO')
  plt.step(bin_centers_lxy, lxy_efficiency1_PAT, where='mid', label='PAT RECO')

  plt.step(bin_centers_lxy, lxy_efficiency2_sct, where='mid', label='Scouting triggers')
  plt.step(bin_centers_lxy, lxy_efficiency2_PAT, where='mid', label='PAT triggers')

  mplhep.cms.label(data=False, com="13.6", loc=0)
  plt.xticks(fontsize=11)
  plt.yticks(fontsize=11)
  plt.xlabel(r"$L_{xy}$ (cm)", fontsize = 16)
  plt.ylabel("Efficiency", fontsize = 16)
  plt.ylim(0, 1.1)
  plt.legend(fontsize=13, loc='top left')
  plt.text(0.35, 0.97, title, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.8))

  plt.text(0.35, 0.93, 'Efficiency of offline reconstruction of 2 gen muons', transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.8))

  plt.grid(True)
  plt.savefig(f"hlt_efficiency_vs_lxy_ctau_{ctau}.png")
  plt.show()
  plt.clf()

  ############################
  #Pt Efficiency
  ############################

  bins_pt = np.linspace(0,40,nbins + 1)
  # Denominator
  pt_counts_sct_denom, _ = np.histogram(pt_sct_with_2genMu, bins=bins_pt)
  pt_counts_PAT_denom, _ = np.histogram(pt_PAT_with_2genMu, bins=bins_pt)
  
  # Numerator
  pt_counts_sct_num1, _ = np.histogram(pt_sct_with_2genMu_and_2sct_mu, bins=bins_pt)
  pt_counts_PAT_num1, _ = np.histogram(pt_PAT_with_2genMu_and_2PAT_mu, bins=bins_pt)

  pt_counts_sct_num2, _ = np.histogram(pt_sct_with_2genMu_and_2sct_mu_and_passTriggers, bins=bins_pt)
  pt_counts_PAT_num2, _ = np.histogram(pt_PAT_with_2genMu_and_2PAT_mu_and_passTriggers, bins=bins_pt)


  pt_efficiency1_sct = pt_counts_sct_num1 / pt_counts_sct_denom
  pt_efficiency1_PAT = pt_counts_PAT_num1 / pt_counts_PAT_denom

  pt_efficiency2_sct = pt_counts_sct_num2 / pt_counts_sct_denom
  pt_efficiency2_PAT = pt_counts_PAT_num2 / pt_counts_PAT_denom


  bin_centers_pt = 0.5 * (bins_pt[:-1] + bins_pt[1:])

  plt.step(bin_centers_pt, pt_efficiency1_sct, where='mid', label='Scouting RECO')
  plt.step(bin_centers_pt, pt_efficiency1_PAT, where='mid', label='PAT RECO')
  plt.step(bin_centers_pt, pt_efficiency2_sct, where='mid', label='Scouting triggers')
  plt.step(bin_centers_pt, pt_efficiency2_PAT, where='mid', label='PAT triggers')

  mplhep.cms.label(data=False, com="13.6", loc=0)
  plt.xticks(fontsize=11)
  plt.yticks(fontsize=11)
  plt.xlabel(r"$P_{t}$ (GeV)", fontsize = 16)
  plt.ylabel("Efficiency", fontsize = 16)
  plt.ylim(0, 1.1)
  plt.legend(fontsize=13, loc='lower center')
  plt.text(0.35, 0.97, title, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.8))

  plt.text(0.35, 0.93, 'Efficiency of offline reconstruction of 2 gen muons', transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.8))

  plt.grid(True)
  plt.savefig(f"hlt_efficiency_vs_pt_ctau_{ctau}.png")
  plt.show()
  plt.clf()
  