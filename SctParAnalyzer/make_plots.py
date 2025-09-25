import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import mplhep
import re
from collections import Counter

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
  
  
  #Scouting vertex variables
  sct_PV_x = tree["sct_PV_x"].array()
  sct_PV_y = tree["sct_PV_y"].array()
  sct_PV_z = tree["sct_PV_z"].array()
  sct_SV_x = tree["sct_SV_x"].array()
  sct_SV_y = tree["sct_SV_y"].array()
  sct_SV_z = tree["sct_SV_z"].array()
  sct_SV_xe = tree["sct_SV_xe"].array()
  sct_SV_ye = tree["sct_SV_ye"].array()
  sct_SV_ze = tree["sct_SV_ze"].array()
  sct_SV_chi2 = tree["sct_SV_chi2"].array()
  sct_SV_prob = tree["sct_SV_prob"].array()
  sct_SV_chi2Ndof = tree["sct_SV_chi2Ndof"].array()
  sct_SV_lxy = tree["sct_SV_lxy"].array()
  sct_SV_l3d = tree["sct_SV_l3d"].array()

  
  #Scouting muon variables:
  sct_Muon_pt_NoVtx = tree["sct_Muon_pt_NoVtx"].array()
  sct_Muon_eta_NoVtx = tree["sct_Muon_eta_NoVtx"].array()
  sct_Muon_phi_NoVtx = tree["sct_Muon_phi_NoVtx"].array()
  nsct_mu_NoVtx = tree["nsct_muons_NoVtx"].array()
  sct_Muon_vtxIdx_NoVtx = tree["sct_Muon_vtxIdx_NoVtx"].array()

  sct_Muon_pt_Vtx = tree["sct_Muon_pt_Vtx"].array()
  sct_Muon_eta_Vtx = tree["sct_Muon_eta_Vtx"].array()
  sct_Muon_phi_Vtx = tree["sct_Muon_phi_Vtx"].array()
  nsct_mu_Vtx = tree["nsct_muons_Vtx"].array()
  sct_Muon_vtxIdx_Vtx = tree["sct_Muon_vtxIdx_Vtx"].array()
  
  ##PAT variables:
  PAT_Muon_pt = tree["PAT_Muon_pt"].array()
  PAT_Muon_eta = tree["PAT_Muon_eta"].array()
  PAT_Muon_phi = tree["PAT_Muon_phi"].array()
  nPAT_mu = tree["nPAT_muons"].array()
  PAT_pair_index = tree["PAT_pair_index"].array()
  PAT_lxy = tree["PAT_lxy"].array()
  PAT_vx = tree["PAT_vx"].array()
  PAT_vy = tree["PAT_vy"].array()

  ##GEN variables:
  Gen_pt = tree["Gen_pt"].array()
  Gen_eta = tree["Gen_eta"].array()
  Gen_phi = tree["Gen_phi"].array()
  #Gen_m = tree["Gen_m"].array()
  Gen_vx = tree["Gen_vx"].array()
  Gen_vy = tree["Gen_vy"].array()
  Gen_vz = tree["Gen_vz"].array()
  Gen_lxy = tree["Gen_lxy"].array()
  #Gen_status = tree["Gen_status"].array()
  #Gen_index = tree["Gen_index"].array()
  #Gen_pdgId = tree["Gen_pdgId"].array()
  #Gen_motherIndex = tree["Gen_motherIndex"].array()
  #Gen_motherPdgId = tree["Gen_motherPdgId"].array()
  #Gen_ct = tree["Gen_ct"].array() 
  ngenMuons = tree["Number_GenMuons"].array()
  '''
  #For debugging purposes
  is_PAT_match = []
  matched_indices_Pat = []
  is_NoVtx_match = []
  matched_indices_NoVtx = []

  is_Vtx_match = []
  matched_indices_Vtx = []
  dr_threshold = 0.1

  for i in range(50):

      event_PAT_match = [False] * len(PAT_Muon_pt[i])
      event_matches_Pat = {}

      event_NoVtx_match = [False] * len(sct_Muon_pt_NoVtx[i])
      event_matches_NoVtx = []

      event_Vtx_match = [False] * len(sct_Muon_pt_Vtx[i])
      event_matches_Vtx = []

      for j in range(len(Gen_pt[i])):
        for k in range(len(PAT_Muon_pt[i])):
          deta_Pat = PAT_Muon_eta[i][k] - Gen_eta[i][j]
          dphi_Pat = PAT_Muon_phi[i][k] - Gen_phi[i][j]
          if dphi_Pat > np.pi:
              dphi_Pat = 2*np.pi - dphi_Pat
          elif dphi_Pat < -np.pi:
              dphi_Pat += 2*np.pi
          dr_Pat = np.sqrt(deta_Pat**2 + dphi_Pat**2)
          if dr_Pat < dr_threshold:
            event_PAT_match[k] = True
            if j not in event_matches_Pat:
                event_matches_Pat[j] = []
            event_matches_Pat[j].append(k)
      print(event_matches_Pat)
      print(len(Gen_pt[i]))
      print(Gen_pt[i])
      print(len(PAT_Muon_pt[i]))
      print(PAT_Muon_pt[i])
      is_PAT_match.append(event_PAT_match)
      matched_indices_Pat.append(event_matches_Pat)
  print(matched_indices_Pat)
  '''
  Draw_plots = True

  if Draw_plots:
    print('Plotting')

    #info_text = f"$m_A$ = {mA} GeV\n$m_{{\pi}}$ = {mpi} GeV\n$c\\tau$ = {ctau} mm"
    info_text = ""

    pt_bins = np.arange(0, 60 + 1, 1)
    nevents = len(Gen_eta)
    
    ################################
    ####### pt distributions #######
    ################################

    selected_gen_pt = []
    selected_Pat_pt = []
    selected_sct_pt_NoVtx = []
    selected_sct_pt_Vtx = []

    subleading_gen_pt = []
    subleading_Pat_pt = []
    subleading_NoVtx_pt = []
    subleading_Vtx_pt = []

    subleading_trigger_Pat_pt = []
    subleading_trigger_NoVtx_pt = []
    subleading_trigger_Vtx_pt = []

    is_PAT_match = []
    matched_indices_Pat = []

    is_NoVtx_match = []
    matched_indices_NoVtx = []

    is_Vtx_match = []
    matched_indices_Vtx = []

    dr_threshold = 0.1

    for i in range(nevents):
      sorted_gen_pts = sorted(Gen_pt[i], reverse=True)
      subleading_gen_pt.append(sorted_gen_pts[1])

      event_PAT_match = [False] * len(PAT_Muon_pt[i])
      event_matches_Pat = {}

      event_NoVtx_match = [False] * len(sct_Muon_pt_NoVtx[i])
      event_matches_NoVtx = {}

      event_Vtx_match = [False] * len(sct_Muon_pt_Vtx[i])
      event_matches_Vtx = {}

      for j in range(len(Gen_pt[i])):
        for k in range(len(PAT_Muon_pt[i])):
          deta_Pat = PAT_Muon_eta[i][k] - Gen_eta[i][j]
          dphi_Pat = PAT_Muon_phi[i][k] - Gen_phi[i][j]
          if dphi_Pat > np.pi:
              dphi_Pat = 2*np.pi - dphi_Pat
          elif dphi_Pat < -np.pi:
              dphi_Pat += 2*np.pi
          dr_Pat = np.sqrt(deta_Pat**2 + dphi_Pat**2)
          if dr_Pat < dr_threshold:
            event_PAT_match[k] = True
            if j not in event_matches_Pat:
                event_matches_Pat[j] = []
            event_matches_Pat[j].append(k)
            


        for k in range(len(sct_Muon_pt_NoVtx[i])):
          deta_NoVtx = sct_Muon_eta_NoVtx[i][k] - Gen_eta[i][j]
          dphi_NoVtx = sct_Muon_phi_NoVtx[i][k] - Gen_phi[i][j]
          if dphi_NoVtx > np.pi:
              dphi_NoVtx = 2*np.pi - dphi_NoVtx
          elif dphi_NoVtx < -np.pi:
              dphi_NoVtx += 2*np.pi
          dr_NoVtx = np.sqrt(deta_NoVtx**2 + dphi_NoVtx**2)
          if dr_NoVtx < dr_threshold:
            event_NoVtx_match[k] = True
            if j not in event_matches_NoVtx:
                event_matches_NoVtx[j] = []
            event_matches_NoVtx[j].append(k) 

        for k in range(len(sct_Muon_pt_Vtx[i])):
          deta_Vtx = sct_Muon_eta_Vtx[i][k] - Gen_eta[i][j]
          dphi_Vtx = sct_Muon_phi_Vtx[i][k] - Gen_phi[i][j]
          if dphi_Vtx > np.pi:
              dphi_Vtx = 2*np.pi - dphi_Vtx
          elif dphi_Vtx < -np.pi:
              dphi_Vtx += 2*np.pi
          dr_Vtx = np.sqrt(deta_Vtx**2 + dphi_Vtx**2)
          if dr_Vtx < dr_threshold:
            event_Vtx_match[k] = True
            if j not in event_matches_Vtx:
                event_matches_Vtx[j] = []
            event_matches_Vtx[j].append(k)
            
      is_PAT_match.append(event_PAT_match)
      matched_indices_Pat.append(event_matches_Pat)

      is_NoVtx_match.append(event_NoVtx_match)
      matched_indices_NoVtx.append(event_matches_NoVtx)

      is_Vtx_match.append(event_Vtx_match)
      matched_indices_Vtx.append(event_matches_Vtx)
      
      ''' 
      matched_pat_pts = [pt for j, pt in enumerate(PAT_Muon_pt[i]) if is_PAT_match[i][j]]
      if len(matched_pat_pts) >= 2:
        subleading_Pat_pt.append(sorted_gen_pts[1])
        if pass_L1_doubleMu[i] and pass_DoubleMu_Parking[i]:
          subleading_trigger_Pat_pt.append(sorted_gen_pts[1])

      matched_NoVtx_pts = [pt for j, pt in enumerate(sct_Muon_pt_NoVtx[i]) if is_NoVtx_match[i][j]]
      if len(matched_NoVtx_pts) >= 2:
        subleading_NoVtx_pt.append(sorted_gen_pts[1])
        if pass_L1_doubleMu[i] and pass_DoubleMu_Scouting[i]:
          subleading_trigger_NoVtx_pt.append(sorted_gen_pts[1])

      matched_Vtx_pts = [pt for j, pt in enumerate(sct_Muon_pt_Vtx[i]) if is_Vtx_match[i][j]]
      if len(matched_Vtx_pts) >= 2:
        subleading_Vtx_pt.append(sorted_gen_pts[1])
        if pass_L1_doubleMu[i] and pass_DoubleMu_Scouting[i]:
          subleading_trigger_Vtx_pt.append(sorted_gen_pts[1])

      for j in range(len(Gen_pt[i])):
        selected_gen_pt.append(Gen_pt[i][j])

      for (gen_idx, pat_idx) in matched_indices_Pat[i]:
        selected_Pat_pt.append(Gen_pt[i][gen_idx])

      for (gen_idx, pat_idx) in matched_indices_NoVtx[i]:
        selected_sct_pt_NoVtx.append(Gen_pt[i][gen_idx])

      for (gen_idx, pat_idx) in matched_indices_Vtx[i]:
          selected_sct_pt_Vtx.append(Gen_pt[i][gen_idx])

    plt.hist(selected_gen_pt, bins=pt_bins, histtype='step', label="Gen pT (≥2 muons)")
    plt.hist(selected_Pat_pt, bins=pt_bins, histtype='step', label="PAT Muon pT")
    plt.hist(selected_sct_pt_NoVtx, bins=pt_bins, histtype='step', label="SCT Muon pT NoVtx")
    plt.hist(selected_sct_pt_Vtx, bins=pt_bins, histtype='step', label="SCT Muon pT Vtx")

    plt.xlabel("Muon pT [GeV]")
    plt.ylabel("Entries")
    plt.legend()
    plt.text(0.75, 0.5, info_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='center', bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.show()
    plt.savefig(f"pt_hist_ctau_{ctau}.png")
    plt.clf()

    ################################
    ######### pt triggers ##########
    ################################
    pt_trigger_bins = np.arange(0, 25 + 1, 1)

    plt.hist(subleading_gen_pt, bins=pt_trigger_bins, histtype='step', label="GEN pT (≥2 muons)")
    plt.hist(subleading_Pat_pt, bins=pt_trigger_bins, histtype='step', label="PAT event RECO")
    plt.hist(subleading_NoVtx_pt, bins=pt_trigger_bins, histtype='step', label="SCT NoVtx event RECO")
    plt.hist(subleading_Vtx_pt, bins=pt_trigger_bins, histtype='step', label="SCT Vtx event RECO")
    plt.hist(subleading_trigger_Pat_pt, bins=pt_trigger_bins, histtype='step', label="PAT passing Double Muon triggers")
    plt.hist(subleading_trigger_NoVtx_pt, bins=pt_trigger_bins, histtype='step', label="SCT NoVtx passing Double Muon triggers")
    plt.hist(subleading_trigger_Vtx_pt, bins=pt_trigger_bins, histtype='step', label="SCT Vtx passing Double Muon triggers")
    plt.xlabel("Subleading Muon pT [GeV]")
    plt.ylabel("Entries")
    plt.legend()
    plt.text(0.75, 0.4, info_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='center', bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.show()
    plt.savefig(f"pt_hist_triggers_ctau_{ctau}.png")
    plt.clf()

    ################################
    ####### pt efficiencies ########
    ################################

    gen_pt_counts, _  = np.histogram(selected_gen_pt, bins=pt_bins)
    pat_pt_counts, _  = np.histogram(selected_Pat_pt, bins=pt_bins)
    sct_pt_counts_NoVtx, _  = np.histogram(selected_sct_pt_NoVtx, bins=pt_bins)
    sct_pt_counts_Vtx, _  = np.histogram(selected_sct_pt_Vtx, bins=pt_bins)

    eff_pt_pat = np.divide(pat_pt_counts, gen_pt_counts, out=np.zeros_like(pat_pt_counts, dtype=float), where=gen_pt_counts!=0)
    eff_pt_sct_NoVtx = np.divide(sct_pt_counts_NoVtx, gen_pt_counts, out=np.zeros_like(sct_pt_counts_NoVtx, dtype=float), where=gen_pt_counts!=0)
    eff_pt_sct_Vtx = np.divide(sct_pt_counts_Vtx, gen_pt_counts, out=np.zeros_like(sct_pt_counts_Vtx, dtype=float), where=gen_pt_counts!=0)

    bin_centers = 0.5 * (pt_bins[1:] + pt_bins[:-1])
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
    plt.clf()

    ################################
    ### pt trigger efficiencies ####
    ################################

    gen_subleading_pt_counts, _ = np.histogram(subleading_gen_pt, bins=pt_trigger_bins)
    pat_event_RECO_pt_counts, _  = np.histogram(subleading_Pat_pt, bins=pt_trigger_bins)
    sct_event_RECO_pt_counts_NoVtx, _  = np.histogram(subleading_NoVtx_pt, bins=pt_trigger_bins)
    sct_event_RECO_pt_counts_Vtx, _  = np.histogram(subleading_Vtx_pt, bins=pt_trigger_bins)

    pat_trigger_pt_counts, _  = np.histogram(subleading_trigger_Pat_pt, bins=pt_trigger_bins)
    sct_trigger_pt_counts_NoVtx, _  = np.histogram(subleading_trigger_NoVtx_pt, bins=pt_trigger_bins)
    sct_trigger_pt_counts_Vtx, _  = np.histogram(subleading_trigger_Vtx_pt, bins=pt_trigger_bins)

    eff_event_RECO_pt_pat = np.divide(pat_event_RECO_pt_counts, gen_subleading_pt_counts, out=np.zeros_like(pat_event_RECO_pt_counts, dtype=float), where=gen_subleading_pt_counts!=0)
    eff_event_RECO_pt_sct_NoVtx = np.divide(sct_event_RECO_pt_counts_NoVtx, gen_subleading_pt_counts, out=np.zeros_like(sct_event_RECO_pt_counts_NoVtx, dtype=float), where=gen_subleading_pt_counts!=0)
    eff_event_RECO_pt_sct_Vtx = np.divide(sct_event_RECO_pt_counts_Vtx, gen_subleading_pt_counts, out=np.zeros_like(sct_event_RECO_pt_counts_Vtx, dtype=float), where=gen_subleading_pt_counts!=0)

    eff_trigger_pt_pat = np.divide(pat_trigger_pt_counts, gen_subleading_pt_counts, out=np.zeros_like(pat_trigger_pt_counts, dtype=float), where=gen_subleading_pt_counts!=0)
    eff_trigger_pt_sct_NoVtx = np.divide(sct_trigger_pt_counts_NoVtx, gen_subleading_pt_counts, out=np.zeros_like(sct_trigger_pt_counts_NoVtx, dtype=float), where=gen_subleading_pt_counts!=0)
    eff_trigger_pt_sct_Vtx = np.divide(sct_trigger_pt_counts_Vtx, gen_subleading_pt_counts, out=np.zeros_like(sct_trigger_pt_counts_Vtx, dtype=float), where=gen_subleading_pt_counts!=0)

    bin_centers = 0.5 * (pt_trigger_bins[1:] + pt_trigger_bins[:-1])
    plt.figure()
    plt.plot(bin_centers, eff_event_RECO_pt_pat, drawstyle='steps-mid', label='PAT Event RECO')
    plt.plot(bin_centers, eff_event_RECO_pt_sct_NoVtx, drawstyle='steps-mid', label='SCT NoVtx Event RECO')
    plt.plot(bin_centers, eff_event_RECO_pt_sct_Vtx, drawstyle='steps-mid', label='SCT Vtx Event RECO')

    plt.plot(bin_centers, eff_trigger_pt_pat, drawstyle='steps-mid', label='Trigger PAT')
    plt.plot(bin_centers, eff_trigger_pt_sct_NoVtx, drawstyle='steps-mid', label='Trigger SCT NoVtx')
    plt.plot(bin_centers, eff_trigger_pt_sct_Vtx, drawstyle='steps-mid', label='Trigger SCT Vtx')

    plt.xlabel("Subleading Muon $p_T$ [GeV]")
    plt.ylabel("Efficiency")
    plt.ylim(0, 1.1)
    plt.grid(True)
    plt.legend()
    plt.savefig(f"pt_trigger_efficiency_ctau_{ctau}.png")
    plt.clf()
    '''

    ################################
    ###### lxy distributions #######
    ################################
    selected_lxy_NoVtx = []
    selected_lxy_Vtx = []
    selected_lxy_Pat = []
    selected_Gen_lxy = []

    selected_trigger_lxy_NoVtx = []
    selected_trigger_lxy_Vtx = []
    selected_trigger_lxy_Pat = []
    selected_trigger_lxy_Gen = []


    for i in range(nevents):
      gen_selected = []
      for j in range(len(Gen_lxy[i])):
        #if Gen_lxy[i][j] not in gen_selected:
          gen_selected.append(Gen_lxy[i][j])
          selected_Gen_lxy.append(Gen_lxy[i][j])
      selected_trigger_lxy_Gen.append(min(gen_selected))

      gen_to_pat_match = matched_indices_Pat[i] if isinstance(matched_indices_Pat[i], dict) else {}
      gen_to_NoVtx_match = matched_indices_NoVtx[i] if isinstance(matched_indices_NoVtx[i], dict) else {}
      gen_to_Vtx_match = matched_indices_Vtx[i] if isinstance(matched_indices_Vtx[i], dict) else {}

      Pat_selected = []
      NoVtx_selected = []
      Vtx_selected = []
      for j1 in range(len(Gen_lxy[i])):
        if j1 in gen_to_pat_match:
          selected_lxy_Pat.append(Gen_lxy[i][j1])
          Pat_selected.append(Gen_lxy[i][j1])
        if j1 in gen_to_NoVtx_match:
          selected_lxy_NoVtx.append(Gen_lxy[i][j1])
          NoVtx_selected.append(Gen_lxy[i][j1])
        if j1 in gen_to_Vtx_match:
          selected_lxy_Vtx.append(Gen_lxy[i][j1])
          Vtx_selected.append(Gen_lxy[i][j1])
      if pass_L1_doubleMu[i] and pass_DoubleMu_Parking[i]:
        if NoVtx_selected:
          selected_trigger_lxy_NoVtx.append(min(NoVtx_selected))
        if Vtx_selected:
          selected_trigger_lxy_Vtx.append(min(Vtx_selected))
        if Pat_selected:
          selected_trigger_lxy_Pat.append(min(Pat_selected))




    lxy_bins = np.linspace(0, 70, 101)  # 100 bins means 101 edges
    plt.figure()
    plt.hist(selected_Gen_lxy, bins=lxy_bins, histtype='step', label="Gen Lxy (≥2 muons)")
    plt.hist(selected_lxy_Pat, bins=lxy_bins, histtype='step', label="RECO PAT Lxy")
    plt.hist(selected_lxy_NoVtx, bins=lxy_bins, histtype='step', label="RECO SCT Lxy (NoVtx)")
    plt.hist(selected_lxy_Vtx, bins=lxy_bins, histtype='step', label="RECO SCT Lxy (Vtx)")
    plt.xlabel("Lxy [cm]")  # or cm depending on your units
    plt.ylabel("Entries")
    plt.legend()
    plt.text(0.75, 0.5, info_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='center', bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.grid(True)
    plt.savefig(f"lxy_hist_ctau_{ctau}.png")
    plt.clf()
    

    ################################
    ###### lxy efficiencies ########
    ################################

    gen_lxy_counts, _ = np.histogram(selected_Gen_lxy, bins=lxy_bins)
    NoVtx_lxy_counts, _ = np.histogram(selected_lxy_NoVtx, bins=lxy_bins)
    Vtx_lxy_counts, _ = np.histogram(selected_lxy_Vtx, bins=lxy_bins)
    Pat_lxy_counts, _ = np.histogram(selected_lxy_Pat, bins=lxy_bins)

    eff_lxy_NoVtx = np.divide(NoVtx_lxy_counts, gen_lxy_counts,out=np.zeros_like(NoVtx_lxy_counts, dtype=float), where=gen_lxy_counts != 0)
    eff_lxy_Vtx = np.divide(Vtx_lxy_counts, gen_lxy_counts,out=np.zeros_like(Vtx_lxy_counts, dtype=float), where=gen_lxy_counts != 0)
    eff_lxy_Pat = np.divide(Pat_lxy_counts, gen_lxy_counts,out=np.zeros_like(Pat_lxy_counts, dtype=float), where=gen_lxy_counts != 0)

    bin_centers = 0.5 * (lxy_bins[:-1] + lxy_bins[1:])

    plt.figure()
    plt.plot(bin_centers, eff_lxy_Pat, drawstyle='steps-mid', label="Pat")
    plt.plot(bin_centers, eff_lxy_NoVtx, drawstyle='steps-mid', label="SCT NoVtx")
    plt.plot(bin_centers, eff_lxy_Vtx, drawstyle='steps-mid', label="SCT Vtx")
    
    plt.xlabel("Lxy [cm]")
    plt.ylabel("Efficiency")
    plt.ylim(0, 1.1)
    plt.grid(True)
    plt.legend()
    plt.savefig(f"lxy_efficiency_ctau_{ctau}.png")
    plt.clf()

    ################################
    ######## lxy triggers ##########
    ################################

    plt.figure()
    plt.hist(selected_trigger_lxy_Gen, bins=lxy_bins, histtype='step', label="Gen Lxy (≥2 muons)")
    plt.hist(selected_trigger_lxy_Pat, bins=lxy_bins, histtype='step', label="PAT Lxy")
    plt.hist(selected_trigger_lxy_NoVtx, bins=lxy_bins, histtype='step', label="SCT Lxy (NoVtx)")
    plt.hist(selected_trigger_lxy_Vtx, bins=lxy_bins, histtype='step', label="SCT Lxy (Vtx)")
    plt.xlabel("Lxy [cm]")
    plt.ylabel("Entries")
    plt.legend()
    plt.text(0.75, 0.5, info_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='center', bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.grid(True)
    plt.savefig(f"lxy_trigger_hist_ctau_{ctau}.png")
    plt.clf()
    

    ################################
    ###### lxy trigger efficiencies ########
    ################################

    gen_trigger_lxy_counts, _ = np.histogram(selected_trigger_lxy_Gen, bins=lxy_bins)
    NoVtx_trigger_lxy_counts, _ = np.histogram(selected_trigger_lxy_NoVtx, bins=lxy_bins)
    Vtx_trigger_lxy_counts, _ = np.histogram(selected_trigger_lxy_Vtx, bins=lxy_bins)
    Pat_trigger_lxy_counts, _ = np.histogram(selected_trigger_lxy_Pat, bins=lxy_bins)

    eff_trigger_lxy_NoVtx = np.divide(NoVtx_trigger_lxy_counts, gen_trigger_lxy_counts,out=np.zeros_like(NoVtx_trigger_lxy_counts, dtype=float), where=gen_trigger_lxy_counts != 0)
    eff_trigger_lxy_Vtx = np.divide(Vtx_trigger_lxy_counts, gen_trigger_lxy_counts,out=np.zeros_like(Vtx_trigger_lxy_counts, dtype=float), where=gen_trigger_lxy_counts != 0)
    eff_trigger_lxy_Pat = np.divide(Pat_trigger_lxy_counts, gen_trigger_lxy_counts,out=np.zeros_like(Pat_trigger_lxy_counts, dtype=float), where=gen_trigger_lxy_counts != 0)


    plt.figure()
    plt.plot(bin_centers, eff_trigger_lxy_Pat, drawstyle='steps-mid', label="Pat")
    plt.plot(bin_centers, eff_trigger_lxy_NoVtx, drawstyle='steps-mid', label="SCT NoVtx")
    plt.plot(bin_centers, eff_trigger_lxy_Vtx, drawstyle='steps-mid', label="SCT NoVtx")
    
    plt.xlabel("Lxy [cm]")
    plt.ylabel("Efficiency")
    plt.ylim(0, 1.1)
    plt.grid(True)
    plt.legend()
    plt.savefig(f"lxy_trigger_efficiency_ctau_{ctau}.png")
    plt.clf()
    
