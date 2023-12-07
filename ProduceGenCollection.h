//Gen Muon
void AnaScript::GenMuon()
{
  int grandmomid=-1;
  for(unsigned int i=0; i< (*nGenPart); i++){
    Particle temp; temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
    temp.status = GenPart_status[i];
    temp.ind = i;
    temp.pdgid = GenPart_pdgId[i];
    temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
    bool passcutmu= abs(temp.pdgid)==13 && temp.status==1 && temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
    bool motherisnotmu= abs(MotherID(i,GenPart_genPartIdxMother[i]))!=13 && GenPart_pdgId[GenPart_genPartIdxMother[i]]!=22;
    passcutmu = passcutmu && motherisnotmu; 
    if(passcutmu){     	
      genMuon.push_back(temp);
      genLightLepton.push_back(temp);
      h.genpltmu[0]->Fill(temp.v.Pt());
      h.genpltmu[1]->Fill(temp.v.Eta());
      h.genpltmu[2]->Fill(temp.v.Phi());
      h.genpltmu[3]->Fill(MotherID(i,GenPart_genPartIdxMother[i]));
    }
  }
}

//Gen Electron
void AnaScript::GenElectron()
{
  int grandmomid=-1;
  for(unsigned int i=0; i< (*nGenPart); i++){
    Particle temp; temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
    temp.status = GenPart_status[i];
    temp.ind = i;
    temp.pdgid = GenPart_pdgId[i];
    temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
    //GenElectron Block
    bool passcutele= abs(temp.pdgid)==11 && temp.status==1 && temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
    bool motherisnotele= abs(MotherID(i,GenPart_genPartIdxMother[i]))!=11 && GenPart_pdgId[GenPart_genPartIdxMother[i]]!=22;
    //Final selection
    passcutele = passcutele && motherisnotele; 
    
    if(passcutele){
      genElectron.push_back(temp);
      genLightLepton.push_back(temp);
      h.genpltele[0]->Fill(temp.v.Pt());
      h.genpltele[1]->Fill(temp.v.Eta());
      h.genpltele[2]->Fill(temp.v.Phi());
      h.genpltele[3]->Fill(MotherID(i,GenPart_genPartIdxMother[i]));
    }
  }
}
//VLL Doublet Yield
void AnaScript::VLLDoublet()
{
  for(unsigned int i=0; i< (*nGenPart); i++){
    Particle temp; temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
    temp.status = GenPart_status[i];
    temp.ind = i;
    temp.pdgid = GenPart_pdgId[i];

    bool passcutVLLDoublet = abs(temp.pdgid)==17;
    bool motherisnotVLLDoublet = abs(GenPart_pdgId[GenPart_genPartIdxMother[i]])!=17; //Mother used directly since it is an unstable particle
    passcutVLLDoublet = passcutVLLDoublet && motherisnotVLLDoublet; 
    if(passcutVLLDoublet)
      genVLLDoublet.push_back(temp);

    bool passcutVLLNuDoublet = abs(temp.pdgid)==18;
    bool motherisnotVLLNuDoublet = abs(GenPart_pdgId[GenPart_genPartIdxMother[i]])!=18; //Mother used directly since it is an unstable particle
    passcutVLLNuDoublet = passcutVLLNuDoublet && motherisnotVLLNuDoublet; 
    if(passcutVLLNuDoublet)
      genVLLNuDoublet.push_back(temp);
  }
}
