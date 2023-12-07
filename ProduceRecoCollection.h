/*************************************************************************
 *                            Light Leptons                              *
 *************************************************************************/

void AnaScript::createLightLeptons(){

  /***********************************
   *              Muons              *
   ***********************************/
  for(unsigned int i=0; i< (*nMuon); i++){
    Particle temp;
    temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105);
    temp.id = -13*Muon_charge[i];
    temp.ind = i;
    temp.charge = Muon_charge[i];
    temp.sip3d = Muon_sip3d[i];
    temp.reliso03 = Muon_pfRelIso03_all[i];
    temp.reliso04 = Muon_pfRelIso04_all[i];

    
    bool ptetacut = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4; 
    bool promptmuon = fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
    //Loose cuts on muons:
    bool passcut_loosemuon = ptetacut && Muon_pfRelIso04_all[i]<0.30 && Muon_looseId[i] && promptmuon;
    //Medium cuts on muons:
    bool passcut_mediummuon = ptetacut && promptmuon && Muon_mediumId[i] && Muon_pfRelIso03_all[i]<0.15;

    //GenMatching Muons
    float dRMin_Muon_GenMuon = 1000.0;
    int index_Mu = -1;
    for(int k=0; k<(int)genMuon.size(); k++){
      float dR_Muon_GenMuon = temp.v.DeltaR(genMuon.at(k).v);
      if(dR_Muon_GenMuon<dRMin_Muon_GenMuon){
	dRMin_Muon_GenMuon = dR_Muon_GenMuon;
	index_Mu = k;
      }
      //Investigating the decay modes for muons in the final state
      int RecoMuon_MotherID = MotherID(genMuon.at(index_Mu).ind,GenPart_genPartIdxMother[genMuon.at(index_Mu).ind]);
      h.momid[0]->Fill(RecoMuon_MotherID);
      if(RecoMuon_MotherID == 23)//Z
	h.decayMode[0]->Fill(1);
      if(RecoMuon_MotherID == 24)//W
	h.decayMode[0]->Fill(2);
      if(RecoMuon_MotherID == 25)//H
	h.decayMode[0]->Fill(3);
      if(RecoMuon_MotherID == 17)//VLL
	h.decayMode[0]->Fill(4);
      if(RecoMuon_MotherID == 18)//VLL_Neutrino
	h.decayMode[0]->Fill(5);
    }
    h.dR[0]->Fill(dRMin_Muon_GenMuon);
    
    if(passcut_mediummuon){
      Muon.push_back(temp);
      LightLepton.push_back(temp);
    }
    if(passcut_loosemuon)
      loosellep.push_back(temp);

  }//for muons

  /***********************************
   *            Electrons            *
   ***********************************/
  for(unsigned int i=0; i< (*nElectron); i++){
    Particle temp;
    temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); 
    temp.id = -11*Electron_charge[i];
    temp.ind = i;
    temp.charge = Electron_charge[i];

    bool ptetacut = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
    bool cleaned_from_muons = clean_from_array(temp, Muon, 0.5);
    bool isprompt = false;
    if(fabs(temp.v.Eta())<=1.479){//for barrel
      if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	isprompt = true;
    }
    if(fabs(temp.v.Eta())>1.479){//endcap
      if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	isprompt = true;
    }
    bool passcut_looseele = ptetacut && Electron_cutBased[i]>1 && isprompt;
    bool passcut_mediumele = ptetacut && isprompt && Electron_cutBased[i]>2 && cleaned_from_muons;

    //GenMatching Electrons
    float dRMin_Electron_GenElectron = 1000.0;
    int index_Ele = -1;
    for(int k=0; k<(int)genElectron.size(); k++){
      float dR_Electron_GenElectron = temp.v.DeltaR(genElectron.at(k).v);
      if(dR_Electron_GenElectron<dRMin_Electron_GenElectron){
	dRMin_Electron_GenElectron = dR_Electron_GenElectron;
	index_Ele = k;
      }
      //Investigating the decay modes for electrons in the final state
      int RecoElectron_MotherID = MotherID(genElectron.at(index_Ele).ind,GenPart_genPartIdxMother[genElectron.at(index_Ele).ind]);
      h.momid[1]->Fill(RecoElectron_MotherID);
      if(RecoElectron_MotherID == 23)//Z
	h.decayMode[1]->Fill(1);
      if(RecoElectron_MotherID == 24)//W
	h.decayMode[1]->Fill(2);
      if(RecoElectron_MotherID == 25)//H
	h.decayMode[1]->Fill(3);
      if(RecoElectron_MotherID == 17)//VLL
	h.decayMode[1]->Fill(4);
      if(RecoElectron_MotherID == 18)//VLL_Neutrino
	h.decayMode[1]->Fill(5);
    }
    h.dR[1]->Fill(dRMin_Electron_GenElectron);
    
    
    if(passcut_mediumele){
      Electron.push_back(temp);
      LightLepton.push_back(temp);
    }
    if(passcut_looseele)
      loosellep.push_back(temp);
  }//For electrons
}

/***********************************
 *              Jets               *
 ***********************************/
void AnaScript::createJets(){
  for(unsigned int i=0; i< (*nJet); i++){
    Particle temp;
    temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
    temp.ind = i;
    
    bool ptetacut = temp.v.Pt()>30 && fabs(temp.v.Eta())<2.4;
    bool cleaned_from_leptons = clean_from_array(temp, loosellep, 0.5);
    bool jetID = _year == 2016 ? Jet_jetId[i]>=1 : Jet_jetId[i]>=2; //if 2016, >=1; else >=2
    bool passcut = ptetacut && cleaned_from_leptons && jetID;
    
    if(passcut){
      Jets.push_back(temp);
      if(Jet_btagDeepB[i]>0.4184)
	bJets.push_back(temp);
    }
  }
}


