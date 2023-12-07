void AnaScript::BookHistograms(){
  h.nevt = new TH1F("nEvents","0-Total events, 1-Total events ran, 2-Total events with trigger applied",5,-1,4);
  h.acceptance = new TH1F("Acceptance","Acceptance_0(N_Generated)_1(4L)_2(3L)_3(2LSS)_4(1L2J)",10,0,10);

  /******************************************************
   *    Object level plots: (before event selection)    *
   ******************************************************/
  //GenMuons
  h.genpltmu[0]  = new TH1F("GenMuon_pt" ,      "GenMuon_pt",  500,0,500);
  h.genpltmu[1]  = new TH1F("GenMuon_eta",      "GenMuon_eta", 200,-4,4);
  h.genpltmu[2]  = new TH1F("GenMuon_phi",      "GenMuon_phi", 200,-4,4);
  h.genpltmu[3]  = new TH1F("GenMuon_MotherID", "GenMuon_MotherID", 1200,-600,600);
  
  //GenElectrons
  h.genpltele[0]  = new TH1F("GenElectron_pt" ,      "GenElectron_pt",  500,0,500);
  h.genpltele[1]  = new TH1F("GenElectron_eta",      "GenElectron_eta", 200,-4,4);
  h.genpltele[2]  = new TH1F("GenElectron_phi",      "GenElectron_phi", 200,-4,4);
  h.genpltele[3]  = new TH1F("GenElectron_MotherID", "GenElectron_MotherID", 1200,-600,600);

  //VLL Doublet
  h.vllD[0]  = new TH1F("nVLLDoublet",   "nVLLDoublet",   5,0,5);
  h.vllNuD[0]= new TH1F("nVLLNuDoublet", "nVLLNuDoublet", 5,0,5);
    
  //Light Leptons
  h.llep[0]  = new TH1F("nLightLeptons", "nLightLeptons", 10,0,10);

  //Electrons
  h.ele[0]  = new TH1F("nElectrons",                  "nElectrons",                  10,0,10);
  h.ele[1]  = new TH1F("Electron_pt",                 "Electron_pt",                 500,0,500);
  h.ele[2]  = new TH1F("Electron_eta",                "Electron_eta",                200,-4,4);
  h.ele[3]  = new TH1F("Electron_phi",                "Electron_phi",                200,-4,4);
  h.ele[4]  = new TH1F("LeadingElectron_pt",          "LeadingElectron_pt",          500,0,500);
  h.ele[5]  = new TH1F("LeadingElectron_eta",         "LeadingElectron_eta",         200,-4,4);
  h.ele[6]  = new TH1F("LeadingElectron_phi",         "LeadingElectron_phi",         200,-4,4);
  h.ele[7]  = new TH1F("SubLeadingElectron_pt",       "SubLeadingElectron_pt",       500,0,500);
  h.ele[8]  = new TH1F("SubLeadingElectron_eta",      "SubLeadingElectron_eta",      200,-4,4);
  h.ele[9]  = new TH1F("SubLeadingElectron_phi",      "SubLeadingElectron_phi",      200,-4,4);
  h.ele[10] = new TH1F("dPhi_LeadingElectron_MET",    "dPhi_LeadingElectron_MET",    200,-4,4);
  h.ele[11] = new TH1F("dPhi_SubLeadingElectron_MET", "dPhi_SubLeadingElectron_MET", 200,-4,4);

  //Muons
  h.mu[0]  = new TH1F("nMuons",                  "nMuons",                  10,0,10);
  h.mu[1]  = new TH1F("Muon_pt",                 "Muon_pt",                 500,0,500);
  h.mu[2]  = new TH1F("Muon_eta",                "Muon_eta",                200,-4,4);
  h.mu[3]  = new TH1F("Muon_phi",                "Muon_phi",                200,-4,4);
  h.mu[4]  = new TH1F("LeadingMuon_pt",          "LeadingMuon_pt",          500,0,500);
  h.mu[5]  = new TH1F("LeadingMuon_eta",         "LeadingMuon_eta",         200,-4,4);
  h.mu[6]  = new TH1F("LeadingMuon_phi",         "LeadingMuon_phi",         200,-4,4);
  h.mu[7]  = new TH1F("SubLeadingMuon_pt",       "SubLeadingMuon_pt",       500,0,500);
  h.mu[8]  = new TH1F("SubLeadingMuon_eta",      "SubLeadingMuon_eta",      200,-4,4);
  h.mu[9]  = new TH1F("SubLeadingMuon_phi",      "SubLeadingMuon_phi",      200,-4,4);
  h.mu[10] = new TH1F("dPhi_LeadingMuon_MET",    "dPhi_LeadingMuon_MET",    200,-4,4);
  h.mu[11] = new TH1F("dPhi_SubLeadingMuon_MET", "dPhi_SubLeadingMuon_MET", 200,-4,4);
  
  //MET
  h.met[0] = new TH1F("MET_pt",  "MET_pt",  500,0,500);
  h.met[1] = new TH1F("MET_phi", "MET_phi", 200,-4,4);

  //EventLevel Variables
  h.lt[0]     = new TH1F("LT",     "LT",     1000,0,1000);
  h.lt_met[0] = new TH1F("LT+MET", "LT+MET", 1000,0,1000);

  //dR Plots
  h.dR[0] = new TH1F("dRMin_Muon_GenMuon",         "dRMin_Muon_GenMuon",         50,0,1);
  h.dR[1] = new TH1F("dRMin_Electron_GenElectron", "dRMin_Electron_GenElectron", 50,0,1);

  //MotherID Plots
  h.momid[0] = new TH1F("MotherId_RecoMuons",     "MotherId_RecoMuons(RecoMuons matched to GenMuons)",             1200,-600,600);
  h.momid[1] = new TH1F("MotherId_RecoElectrons", "MotherId_RecoElectrons(RecoElectrons matched to GenElectrons)", 1200,-600,600);

  //Decay Modes
  h.decayMode[0] = new TH1F("RecoMuon_DecayMode",        "RecoMuon_DecayModes_Z(1)_W(2)_H(3)_L(4)_N(5)",     10,0,10);
  h.decayMode[1] = new TH1F("RecoElectron_DecayMode",    "RecoElectron_DecayModes_Z(1)_W(2)_H(3)_L(4)_N(5)", 10,0,10);
  h.decayMode[2] = new TH1F("VLLDoublet_ProductionMode", "VLLDoublet_ProductionModes_LL(1)_LN(2)_NN(3)",     5,0,5);

  //4L Event Variables
  h.hist_l4[0]  = new TH1F("4L_nLightLeptons", "4L_nLightLeptons", 5,0,5);
  h.hist_l4[1]  = new TH1F("4L_L0_pt",         "4L_L0_pt",         500,0,500);
  h.hist_l4[2]  = new TH1F("4L_L1_pt",         "4L_L1_pt",         500,0,500);
  h.hist_l4[3]  = new TH1F("4L_L2_pt",         "4L_L2_pt",         500,0,500);
  h.hist_l4[4]  = new TH1F("4L_L3_pt",         "4L_L3_pt",         500,0,500);
  h.hist_l4[5]  = new TH1F("4L_L0_eta",        "4L_L0_eta",        200,-4,4);
  h.hist_l4[6]  = new TH1F("4L_L1_eta",        "4L_L1_eta",        200,-4,4);
  h.hist_l4[7]  = new TH1F("4L_L2_eta",        "4L_L2_eta",        200,-4,4);
  h.hist_l4[8]  = new TH1F("4L_L3_eta",        "4L_L3_eta",        200,-4,4);
  h.hist_l4[9]  = new TH1F("4L_L0_phi",        "4L_L0_phi",        200,-4,4);
  h.hist_l4[10] = new TH1F("4L_L1_phi",        "4L_L1_phi",        200,-4,4);
  h.hist_l4[11] = new TH1F("4L_L2_phi",        "4L_L2_phi",        200,-4,4);
  h.hist_l4[12] = new TH1F("4L_L3_phi",        "4L_L3_phi",        200,-4,4);
  h.hist_l4[13] = new TH1F("4L_LT",            "4L_LT",            500,0,500);
  h.hist_l4[14] = new TH1F("4L_MET",           "4L_MET",           500,0,500);
  h.hist_l4[15] = new TH1F("4L_LT+MET",        "4L_LT+MET",        500,0,500);
  h.hist_l4[16] = new TH1F("4L_dPhi_L0MET",    "4L_dPhi_L0MET",    100,0,4);
  h.hist_l4[17] = new TH1F("4L_dPhi_L1MET",    "4L_dPhi_L1MET",    100,0,4);
  h.hist_l4[18] = new TH1F("4L_dPhi_L2MET",    "4L_dPhi_L2MET",    100,0,4);
  h.hist_l4[19] = new TH1F("4L_dPhi_L3MET",    "4L_dPhi_L3MET",    100,0,4);
  h.hist_l4[20] = new TH1F("4L_dPhi_L0L1",     "4L_dPhi_L0L1",     100,0,4);
  h.hist_l4[21] = new TH1F("4L_dPhi_L0L2",     "4L_dPhi_L0L2",     100,0,4);
  h.hist_l4[22] = new TH1F("4L_dPhi_L0L3",     "4L_dPhi_L0L3",     100,0,4);
  h.hist_l4[23] = new TH1F("4L_dPhi_L1L2",     "4L_dPhi_L1L2",     100,0,4);
  h.hist_l4[24] = new TH1F("4L_dPhi_L1L3",     "4L_dPhi_L1L3",     100,0,4);
  h.hist_l4[25] = new TH1F("4L_dPhi_L2L3",     "4L_dPhi_L2L3",     100,0,4);
  h.hist_l4[26] = new TH1F("4L_MT_L0",         "4L_MT_L0",         200,0,200);
  h.hist_l4[27] = new TH1F("4L_MT_L1",         "4L_MT_L1",         200,0,200);
  h.hist_l4[28] = new TH1F("4L_MT_L2",         "4L_MT_L2",         200,0,200);
  h.hist_l4[29] = new TH1F("4L_MT_L3",         "4L_MT_L3",         200,0,200);
  h.hist_l4[30] = new TH1F("4L_dR_L0L1",       "4L_dR_L0L1",       100,0,4);
  h.hist_l4[31] = new TH1F("4L_dR_L0L2",       "4L_dR_L0L2",       100,0,4);
  h.hist_l4[32] = new TH1F("4L_dR_L0L3",       "4L_dR_L0L3",       100,0,4);
  h.hist_l4[33] = new TH1F("4L_dR_L1L2",       "4L_dR_L1L2",       100,0,4);
  h.hist_l4[34] = new TH1F("4L_dR_L1L3",       "4L_dR_L1L3",       100,0,4);
  h.hist_l4[35] = new TH1F("4L_dR_L2L3",       "4L_dR_L2L3",       100,0,4);
  h.hist_l4[36] = new TH1F("4L_InvMass_L0L1",  "4L_InvMass_L0L1",  200,0,200);
  h.hist_l4[37] = new TH1F("4L_InvMass_L0L2",  "4L_InvMass_L0L2",  200,0,200);
  h.hist_l4[38] = new TH1F("4L_InvMass_L0L3",  "4L_InvMass_L0L3",  200,0,200);
  h.hist_l4[39] = new TH1F("4L_InvMass_L1L2",  "4L_InvMass_L1L2",  200,0,200);
  h.hist_l4[40] = new TH1F("4L_InvMass_L1L3",  "4L_InvMass_L1L3",  200,0,200);
  h.hist_l4[41] = new TH1F("4L_InvMass_L2L3",  "4L_InvMass_L2L3",  200,0,200);
    

  //3L Event Variables
  h.hist_l3[0]  = new TH1F("3L_nLightLeptons", "3L_nLightLeptons", 5,0,5);
  h.hist_l3[1]  = new TH1F("3L_L0_pt",         "3L_L0_pt",         500,0,500);
  h.hist_l3[2]  = new TH1F("3L_L1_pt",         "3L_L1_pt",         500,0,500);
  h.hist_l3[3]  = new TH1F("3L_L2_pt",         "3L_L2_pt",         500,0,500);
  h.hist_l3[4]  = new TH1F("3L_L0_eta",        "3L_L0_eta",        200,-4,4);
  h.hist_l3[5]  = new TH1F("3L_L1_eta",        "3L_L1_eta",        200,-4,4);
  h.hist_l3[6]  = new TH1F("3L_L2_eta",        "3L_L2_eta",        200,-4,4);
  h.hist_l3[7]  = new TH1F("3L_L0_phi",        "3L_L0_phi",        200,-4,4);
  h.hist_l3[8]  = new TH1F("3L_L1_phi",        "3L_L1_phi",        200,-4,4);
  h.hist_l3[9]  = new TH1F("3L_L2_phi",        "3L_L2_phi",        200,-4,4);
  h.hist_l3[10] = new TH1F("3L_LT",            "3L_LT",            500,0,500);
  h.hist_l3[11] = new TH1F("3L_MET",           "3L_MET",           500,0,500);
  h.hist_l3[12] = new TH1F("3L_LT+MET",        "3L_LT+MET",        500,0,500);
  h.hist_l3[13] = new TH1F("3L_dPhi_L0MET",    "3L_dPhi_L0MET",    100,0,4);
  h.hist_l3[14] = new TH1F("3L_dPhi_L1MET",    "3L_dPhi_L1MET",    100,0,4);
  h.hist_l3[15] = new TH1F("3L_dPhi_L2MET",    "3L_dPhi_L2MET",    100,0,4);
  h.hist_l3[16] = new TH1F("3L_dPhi_L0L1",     "3L_dPhi_L0L1",     100,0,4);
  h.hist_l3[17] = new TH1F("3L_dPhi_L0L2",     "3L_dPhi_L0L2",     100,0,4);
  h.hist_l3[18] = new TH1F("3L_dPhi_L1L2",     "3L_dPhi_L1L2",     100,0,4);
  h.hist_l3[19] = new TH1F("3L_MT_L0",         "3L_MT_L0",         200,0,200);
  h.hist_l3[20] = new TH1F("3L_MT_L1",         "3L_MT_L1",         200,0,200);
  h.hist_l3[21] = new TH1F("3L_MT_L2",         "3L_MT_L2",         200,0,200);
  h.hist_l3[22] = new TH1F("3L_dR_L0L1",       "3L_dR_L0L1",       100,0,4);
  h.hist_l3[23] = new TH1F("3L_dR_L0L2",       "3L_dR_L0L2",       100,0,4);
  h.hist_l3[24] = new TH1F("3L_dR_L1L2",       "3L_dR_L1L2",       100,0,4);
  h.hist_l3[25] = new TH1F("3L_InvMass_L0L1",  "3L_InvMass_L0L1",  200,0,200);
  h.hist_l3[26] = new TH1F("3L_InvMass_L0L2",  "3L_InvMass_L0L2",  200,0,200);
  h.hist_l3[27] = new TH1F("3L_InvMass_L1L2",  "3L_InvMass_L1L2",  200,0,200);

  //2L SS Event Variables
  h.hist_l2SS[0]  = new TH1F("2LSS_nLightLeptons", "2LSS_nLightLeptons", 5,0,5);
  h.hist_l2SS[1]  = new TH1F("2LSS_L0_pt",         "2LSS_L0_pt",         500,0,500);
  h.hist_l2SS[2]  = new TH1F("2LSS_L1_pt",         "2LSS_L1_pt",         500,0,500);
  h.hist_l2SS[3]  = new TH1F("2LSS_L0_eta",        "2LSS_L0_eta",        200,-4,4);
  h.hist_l2SS[4]  = new TH1F("2LSS_L1_eta",        "2LSS_L1_eta",        200,-4,4);
  h.hist_l2SS[5]  = new TH1F("2LSS_L0_phi",        "2LSS_L0_phi",        200,-4,4);
  h.hist_l2SS[6]  = new TH1F("2LSS_L1_phi",        "2LSS_L1_phi",        200,-4,4);
  h.hist_l2SS[7]  = new TH1F("2LSS_LT",            "2LSS_LT",            500,0,500);
  h.hist_l2SS[8]  = new TH1F("2LSS_MET",           "2LSS_MET",           500,0,500);
  h.hist_l2SS[9]  = new TH1F("2LSS_LT+MET",        "2LSS_LT+MET",        500,0,500);
  h.hist_l2SS[10] = new TH1F("2LSS_dPhi_L0MET",    "2LSS_dPhi_L0MET",    100,0,4);
  h.hist_l2SS[11] = new TH1F("2LSS_dPhi_L1MET",    "2LSS_dPhi_L1MET",    100,0,4);
  h.hist_l2SS[12] = new TH1F("2LSS_dPhi_L0L1",     "2LSS_dPhi_L0L1",     100,0,4);
  h.hist_l2SS[13] = new TH1F("2LSS_MT_L0",         "2LSS_MT_L0",         200,0,200);
  h.hist_l2SS[14] = new TH1F("2LSS_MT_L1",         "2LSS_MT_L1",         200,0,200);
  h.hist_l2SS[15] = new TH1F("2LSS_dR_L0L1",       "2LSS_dR_L0L1",       100,0,4);
  h.hist_l2SS[16] = new TH1F("2LSS_InvMass_L0L1",  "2LSS_InvMass_L0L1",  200,0,200);

  //1L2J Event Variables
  h.hist_l12j[0]  = new TH1F("1L2J_nLightLeptons", "1L2J_nLightLeptons", 5,0,5);
  h.hist_l12j[1]  = new TH1F("1L2J_nJets",         "1L2J_nJets",         5,0,5);
  h.hist_l12j[2]  = new TH1F("1L2J_L0_pt",         "1L2J_L0_pt",         500,0,500);
  h.hist_l12j[3]  = new TH1F("1L2J_J0_pt",         "1L2J_J0_pt",         500,0,500);
  h.hist_l12j[4]  = new TH1F("1L2J_J1_pt",         "1L2J_J1_pt",         500,0,500);
  h.hist_l12j[5]  = new TH1F("1L2J_L0_eta",        "1L2J_L0_eta",        200,-4,4);
  h.hist_l12j[6]  = new TH1F("1L2J_J0_eta",        "1L2J_J0_eta",        200,-4,4);
  h.hist_l12j[7]  = new TH1F("1L2J_J1_eta",        "1L2J_J1_eta",        200,-4,4);
  h.hist_l12j[8]  = new TH1F("1L2J_L0_phi",        "1L2J_L0_phi",        200,-4,4);
  h.hist_l12j[9]  = new TH1F("1L2J_J0_phi",        "1L2J_J0_phi",        200,-4,4);
  h.hist_l12j[10] = new TH1F("1L2J_J1_phi",        "1L2J_J1_phi",        200,-4,4);
  h.hist_l12j[11] = new TH1F("1L2J_MT_L0",         "1L2J_MT_L0",         200,0,200);
  h.hist_l12j[12] = new TH1F("1L2J_MT_J0",         "1L2J_MT_J0",         200,0,200);
  h.hist_l12j[13] = new TH1F("1L2J_MT_J1",         "1L2J_MT_J1",         200,0,200);
  h.hist_l12j[14] = new TH1F("1L2J_dPhi_L0MET",    "1L2J_dPhi_L0MET",    100,0,4);
  h.hist_l12j[15] = new TH1F("1L2J_dPhi_J0MET",    "1L2J_dPhi_J0MET",    100,0,4);
  h.hist_l12j[16] = new TH1F("1L2J_dPhi_J1MET",    "1L2J_dPhi_J1MET",    100,0,4);
  h.hist_l12j[17] = new TH1F("1L2J_dPhi_L0J0",     "1L2J_dPhi_L0J0",     100,0,4);
  h.hist_l12j[18] = new TH1F("1L2J_dPhi_L0J1",     "1L2J_dPhi_L0J1",     100,0,4);
  h.hist_l12j[19] = new TH1F("1L2J_dPhi_J0J1",     "1L2J_dPhi_J0J1",     100,0,4);
  h.hist_l12j[20] = new TH1F("1L2J_dPhi_JOJ1_L0",  "1L2J_dPhi_J0J1_L0",  100,0,4);
  h.hist_l12j[21] = new TH1F("1L2J_dPhi_J0J1_MET", "1L2J_dPhi_J0J1_MET", 100,0,4);
  h.hist_l12j[22] = new TH1F("1L2J_dPhi_LJJ_MET",  "1L2J_dPhi_LJJ_MET",  100,0,4);
  h.hist_l12j[23] = new TH1F("1L2J_dPhi_LJJ_L0",   "1L2J_dPhi_LJJ_L0",   100,0,4);
  h.hist_l12j[24] = new TH1F("1L2J_dPhi_LJJ_J0",   "1L2J_dPhi_LJJ_J0",   100,0,4);
  h.hist_l12j[25] = new TH1F("1L2J_dPhi_LJJ_J1",   "1L2J_dPhi_LJJ_J1",   100,0,4);
  h.hist_l12j[26] = new TH1F("1L2J_InvMass_J0J1",  "1L2J_InvMass_J0J1",  200,0,200);
  h.hist_l12j[27] = new TH1F("1L2J_pT_J0J1",       "1L2J_pT_J0J1",       500,0,500);
  h.hist_l12j[28] = new TH1F("1L2J_MT_J0J1",       "1L2J_MT_J0J1",       200,0,200);
  h.hist_l12j[29] = new TH1F("1L2J_dR_J0J1",       "1L2J_dR_J0J1",       100,0,4);
  h.hist_l12j[30] = new TH1F("1L2J_InvMass_LJJ",   "1L2J_InvMass_LJJ",   500,0,500);
  h.hist_l12j[31] = new TH1F("1L2J_pT_LJJ",        "1L2J_pT_LJJ",        500,0,500);
  h.hist_l12j[32] = new TH1F("1L2J_MET",           "1L2J_MET",           500,0,500);
  h.hist_l12j[33] = new TH1F("1L2J_HT",            "1L2J_HT",            500,0,500);
  h.hist_l12j[34] = new TH1F("1L2J_ST",            "1L2J_ST",            1000,0,1000);
  
  //WZ Control Region Plots
  h.WZ_CR[0]  = new TH1F("3M_Chosen12_Mass01",      "3M_Chosen12_Mass01",      200,0,200);
  h.WZ_CR[1]  = new TH1F("3M_Chosen12_TransMass2",  "3M_Chosen12_TransMass2",  200,0,200);
  h.WZ_CR[2]  = new TH1F("3M_Chosen12_Mass02",      "3M_Chosen12_Mass02",      200,0,200);
  h.WZ_CR[3]  = new TH1F("3M_Chosen12_TransMass1",  "3M_Chosen12_TransMass1",  200,0,200);
  h.WZ_CR[4]  = new TH1F("3M_Chosen23_Mass02",      "3M_Chosen23_Mass02",      200,0,200);
  h.WZ_CR[5]  = new TH1F("3M_Chosen23_TransMass1",  "3M_Chosen23_TransMass1",  200,0,200);
  h.WZ_CR[6]  = new TH1F("3M_Chosen23_Mass12",      "3M_Chosen23_Mass12",      200,0,200);
  h.WZ_CR[7]  = new TH1F("3M_Chosen23_TransMass0",  "3M_Chosen23_TransMass0",  200,0,200);
  h.WZ_CR[8]  = new TH1F("3M_Chosen13_Mass01",      "3M_Chosen13_Mass01",      200,0,200);
  h.WZ_CR[9]  = new TH1F("3M_Chosen13_TransMass2",  "3M_Chosen13_TransMass2",  200,0,200);
  h.WZ_CR[10] = new TH1F("3M_Chosen13_Mass12",      "3M_Chosen13_Mass12",      200,0,200);
  h.WZ_CR[11] = new TH1F("3M_Chosen13_TransMass0",  "3M_Chosen13_TransMass0",  200,0,200);
  h.WZ_CR[12] = new TH1F("3E_Chosen12_Mass01",      "3E_Chosen12_Mass01",      200,0,200);
  h.WZ_CR[13] = new TH1F("3E_Chosen12_TransMass2",  "3E_Chosen12_TransMass2",  200,0,200);
  h.WZ_CR[14] = new TH1F("3E_Chosen12_Mass02",      "3E_Chosen12_Mass02",      200,0,200);
  h.WZ_CR[15] = new TH1F("3E_Chosen12_TransMass1",  "3E_Chosen12_TransMass1",  200,0,200);
  h.WZ_CR[16] = new TH1F("3E_Chosen23_Mass02",      "3E_Chosen23_Mass02",      200,0,200);
  h.WZ_CR[17] = new TH1F("3E_Chosen23_TransMass1",  "3E_Chosen23_TransMass1",  200,0,200);
  h.WZ_CR[18] = new TH1F("3E_Chosen23_Mass12",      "3E_Chosen23_Mass12",      200,0,200);
  h.WZ_CR[19] = new TH1F("3E_Chosen23_TransMass0",  "3E_Chosen23_TransMass0",  200,0,200);
  h.WZ_CR[20] = new TH1F("3E_Chosen13_Mass01",      "3E_Chosen13_Mass01",      200,0,200);
  h.WZ_CR[21] = new TH1F("3E_Chosen13_TransMass2",  "3E_Chosen13_TransMass2",  200,0,200);
  h.WZ_CR[22] = new TH1F("3E_Chosen13_Mass12",      "3E_Chosen13_Mass12",      200,0,200);
  h.WZ_CR[23] = new TH1F("3E_Chosen13_TransMass0",  "3E_Chosen13_TransMass0",  200,0,200);
  h.WZ_CR[24] = new TH1F("1E2M_Chosen1_Mass01",     "1E2M_Chosen1_Mass01",     200,0,200);
  h.WZ_CR[25] = new TH1F("1E2M_Chosen1_TransMass2", "1E2M_Chosen1_TransMass2", 200,0,200);
  h.WZ_CR[26] = new TH1F("1E2M_Chosen2_Mass02",     "1E2M_Chosen2_Mass02",     200,0,200);
  h.WZ_CR[27] = new TH1F("1E2M_Chosen2_TransMass1", "1E2M_Chosen2_TransMass1", 200,0,200);
  h.WZ_CR[28] = new TH1F("1E2M_Chosen3_Mass12",     "1E2M_Chosen3_Mass12",     200,0,200);
  h.WZ_CR[29] = new TH1F("1E2M_Chosen3_TransMass0", "1E2M_Chosen3_TransMass0", 200,0,200);
  h.WZ_CR[30] = new TH1F("1M2E_Chosen1_Mass01",     "1M2E_Chosen1_Mass01",     200,0,200);
  h.WZ_CR[31] = new TH1F("1M2E_Chosen1_TransMass2", "1M2E_Chosen1_TransMass2", 200,0,200);
  h.WZ_CR[32] = new TH1F("1M2E_Chosen2_Mass02",     "1M2E_Chosen2_Mass02",     200,0,200);
  h.WZ_CR[33] = new TH1F("1M2E_Chosen2_TransMass1", "1M2E_Chosen2_TransMass1", 200,0,200);
  h.WZ_CR[34] = new TH1F("1M2E_Chosen3_Mass12",     "1M2E_Chosen3_Mass12",     200,0,200);
  h.WZ_CR[35] = new TH1F("1M2E_Chosen3_TransMass0", "1M2E_Chosen3_TransMass0", 200,0,200);
  h.WZ_CR[36] = new TH1F("MassArray_3L_Size",       "MassArray_3L_Size",       5,0,5);
  h.WZ_CR[37] = new TH1F("TrMassArray_3L_Size",     "TrMassArray_3L_Size",     5,0,5);
  h.WZ_CR[38] = new TH1F("WZ_CR_Test",              "WZ_CR_Test",              3,0,3);
}
