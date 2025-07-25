import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import mplhep
import re

#filenames = ["output_ctau-1-mA-2p00-mpi-10.root"]
#filenames = ["output_ctau-10-mA-2p00-mpi-10.root"]
filenames = ["output_ctau-100-mA-2p00-mpi-10.root"]
#filenames = ["output_ctau-1-mA-2p00-mpi-10.root", "output_ctau-10-mA-2p00-mpi-10.root", "output_ctau-100-mA-2p00-mpi-10.root"]
#filenames = ["output_test.root"]
for filename in filenames:
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
  '''
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
  '''
  
  #Scouting vertex variables
  #sct_PV_x = tree["sct_PV_x"].array()
  #sct_PV_y = tree["sct_PV_y"].array()
  #sct_PV_z = tree["sct_PV_z"].array()
  #sct_SV_x = tree["sct_SV_x"].array()
  #sct_SV_y = tree["sct_SV_y"].array()
  #sct_SV_z = tree["sct_SV_z"].array()
  #sct_SV_xe = tree["sct_SV_xe"].array()
  #sct_SV_ye = tree["sct_SV_ye"].array()
  #sct_SV_ze = tree["sct_SV_ze"].array()
  #sct_SV_chi2 = tree["sct_SV_chi2"].array()
  #sct_SV_prob = tree["sct_SV_prob"].array()
  #sct_SV_chi2Ndof = tree["sct_SV_chi2Ndof"].array()
  #sct_SV_lxy = tree["sct_SV_lxy"].array()
  #sct_SV_l3d = tree["sct_SV_l3d"].array()

  
  #Scouting muon variables:
  
  sct_Muon_pt_NoVtx = tree["sct_Muon_pt_NoVtx"].array()
  sct_Muon_eta_NoVtx = tree["sct_Muon_eta_NoVtx"].array()
  sct_Muon_phi_NoVtx = tree["sct_Muon_phi_NoVtx"].array()
  nsct_mu_NoVtx = tree["nsct_muons_NoVtx"].array()

  sct_Muon_pt_Vtx = tree["sct_Muon_pt_Vtx"].array()
  sct_Muon_eta_Vtx = tree["sct_Muon_eta_Vtx"].array()
  sct_Muon_phi_Vtx = tree["sct_Muon_phi_Vtx"].array()
  nsct_mu_Vtx = tree["nsct_muons_Vtx"].array()
  

  
  ##PAT variables:
  PAT_Muon_pt = tree["PAT_Muon_pt"].array()
  PAT_Muon_eta = tree["PAT_Muon_eta"].array()
  PAT_Muon_phi = tree["PAT_Muon_phi"].array()
  nPAT_mu = tree["nPAT_muons"].array()

  ##GEN variables:
  Gen_pt = tree["Gen_pt"].array()
  Gen_eta = tree["Gen_eta"].array()
  Gen_phi = tree["Gen_phi"].array()
  #Gen_m = tree["Gen_m"].array()
  #Gen_vx = tree["Gen_vx"].array()
  #Gen_vy = tree["Gen_vy"].array()
  #Gen_vz = tree["Gen_vz"].array()
  #Gen_lxy = tree["Gen_lxy"].array()
  #Gen_status = tree["Gen_status"].array()
  #Gen_index = tree["Gen_index"].array()
  #Gen_pdgId = tree["Gen_pdgId"].array()
  #Gen_motherIndex = tree["Gen_motherIndex"].array()
  #Gen_motherPdgId = tree["Gen_motherPdgId"].array()
  #Gen_ct = tree["Gen_ct"].array()
  
  ngenMuons = tree["Number_GenMuons"].array()


  Draw_plots = True
  if Draw_plots:
    print('Plotting')

    gen_eta = ak.Array(Gen_eta)
    gen_phi = ak.Array(Gen_phi)

    pat_eta = ak.Array(PAT_Muon_eta)
    pat_phi = ak.Array(PAT_Muon_phi)
    sct_eta_NoVtx = ak.Array(sct_Muon_eta_NoVtx)
    sct_phi_NoVtx = ak.Array(sct_Muon_phi_NoVtx)
    sct_eta_Vtx = ak.Array(sct_Muon_eta_Vtx)
    sct_phi_Vtx = ak.Array(sct_Muon_phi_Vtx)

    deta_Pat = pat_eta[:, :, None] - gen_eta[:, None, :]
    dphi_Pat = pat_phi[:, :, None] - gen_phi[:, None, :]
    deta_sct_NoVtx = sct_eta_NoVtx[:, :, None] - gen_eta[:, None, :]
    dphi_sct_NoVtx = sct_phi_NoVtx[:, :, None] - gen_phi[:, None, :]
    deta_sct_Vtx = sct_eta_Vtx[:, :, None] - gen_eta[:, None, :]
    dphi_sct_Vtx = sct_phi_Vtx[:, :, None] - gen_phi[:, None, :]
    #dphi = (dphi + np.pi) % (2 * np.pi) - np.pi

    dr_Pat = np.sqrt(deta_Pat**2 + dphi_Pat**2)
    dr_sct_NoVtx = np.sqrt(deta_sct_NoVtx**2 + dphi_sct_NoVtx**2)
    dr_sct_Vtx = np.sqrt(deta_sct_Vtx**2 + dphi_sct_Vtx**2)

    dr_threshold = 0.1

    is_PAT_match = ak.any(dr_Pat < dr_threshold, axis=2) 
    is_sct_NoVtx_match = ak.any(dr_sct_NoVtx < dr_threshold, axis=2) 
    is_sct_Vtx_match = ak.any(dr_sct_Vtx < dr_threshold, axis=2) 

    x_bins = np.arange(0, 100 + 1, 1)
    nevents = len(Gen_eta)

    ################################
    ####### pt distributions #######
    ################################

    selected_gen_pt = []
    selected_Pat_pt = []
    selected_sct_pt_NoVtx = []
    selected_sct_pt_Vtx = []
    for i in range(nevents):
      if ngenMuons[i]>=2:
        for j in range(len(Gen_pt[i])):
          selected_gen_pt.append(Gen_pt[i][j])
        for j in range(len(PAT_Muon_pt[i])):
          if is_PAT_match[i][j]:
            selected_Pat_pt.append(PAT_Muon_pt[i][j])
        for j in range(len(sct_Muon_pt_NoVtx[i])):
          if is_sct_NoVtx_match[i][j]:
            selected_sct_pt_NoVtx.append(sct_Muon_pt_NoVtx[i][j])
        for j in range(len(sct_Muon_pt_Vtx[i])):
          if is_sct_Vtx_match[i][j]:
            selected_sct_pt_Vtx.append(sct_Muon_pt_Vtx[i][j])
    
    
    plt.hist(selected_gen_pt, bins=x_bins, histtype='step', label="Gen pT (≥2 muons)")
    plt.hist(selected_Pat_pt, bins=x_bins, histtype='step', label="Matched PAT Muon pT")
    plt.hist(selected_sct_pt_NoVtx, bins=x_bins, histtype='step', label="Matched SCT Muon pT NoVtx")
    plt.hist(selected_sct_pt_Vtx, bins=x_bins, histtype='step', label="Matched SCT Muon pT Vtx")
    plt.xlabel("Muon pT [GeV]")
    plt.ylabel("Entries")
    plt.legend()
    plt.show()
    plt.savefig(f"pt_hist__ctau_{ctau}.png")

    ################################
    ####### pt efficiencies ########
    ################################

    gen_pt_counts, _  = np.histogram(selected_gen_pt, bins=x_bins)
    pat_pt_counts, _  = np.histogram(selected_Pat_pt, bins=x_bins)
    sct_pt_counts_NoVtx, _  = np.histogram(selected_sct_pt_NoVtx, bins=x_bins)
    sct_pt_counts_Vtx, _  = np.histogram(selected_sct_pt_Vtx, bins=x_bins)

    eff_pt_pat = np.divide(pat_pt_counts, gen_pt_counts, out=np.zeros_like(pat_pt_counts, dtype=float), where=gen_pt_counts!=0)
    eff_pt_sct_NoVtx = np.divide(sct_pt_counts_NoVtx, gen_pt_counts, out=np.zeros_like(sct_pt_counts_NoVtx, dtype=float), where=gen_pt_counts!=0)
    eff_pt_sct_Vtx = np.divide(sct_pt_counts_Vtx, gen_pt_counts, out=np.zeros_like(sct_pt_counts_Vtx, dtype=float), where=gen_pt_counts!=0)

    bin_centers = 0.5 * (x_bins[1:] + x_bins[:-1])
    plt.figure()
    plt.plot(bin_centers, eff_pt_pat, drawstyle='steps-mid', label='PAT Efficiency')
    plt.plot(bin_centers, eff_pt_sct_NoVtx, drawstyle='steps-mid', label='SCT NoVtx Efficiency')
    plt.plot(bin_centers, eff_pt_sct_Vtx, drawstyle='steps-mid', label='SCT Vtx Efficiency')
    plt.xlabel("Muon $p_T$ [GeV]")
    plt.ylabel("Efficiency")
    plt.ylim(0, 1.1)
    plt.grid(True)
    plt.legend()
    plt.savefig(f"pt_efficiency_ctau_{ctau}.png")







         

    
  