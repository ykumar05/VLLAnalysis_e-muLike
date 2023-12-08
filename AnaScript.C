#define AnaScript_cxx
#include "AnaScript.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

//Including the header files:
#include "/home/work/ykumar1/Work/VLLAnalysis_e-muLike/AnalyzerCode/BookHistograms.h"
#include "/home/work/ykumar1/Work/VLLAnalysis_e-muLike/AnalyzerCode/CustomFunctions.h"
#include "/home/work/ykumar1/Work/VLLAnalysis_e-muLike/AnalyzerCode/ProduceGenCollection.h"
#include "/home/work/ykumar1/Work/VLLAnalysis_e-muLike/AnalyzerCode/ProduceRecoCollection.h"
#include "/home/work/ykumar1/Work/VLLAnalysis_e-muLike/AnalyzerCode/CR_WZ.h"

void AnaScript::Begin(TTree * /*tree*/)
{
  TString option = GetOption();
}

void AnaScript::SlaveBegin(TTree * /*tree*/)
{
  time(&start);

  TString option = GetOption();
  nEvtTotal = 0;
  nEvtRan = 0;
  nEvtTrigger=0;
  n_l2OS = 0;
  n_l2SS = 0;
  n_l3 = 0;
  n_l4 = 0;
  n_l12j = 0;

  _HstFile = new TFile(_HstFileName,"recreate");
  //Call the function to book the histograms we declared in Hists.
  BookHistograms();
}

void AnaScript::SlaveTerminate()
{
  _HstFile->Write();
  _HstFile->Close();

  float goodevtfrac = ((float)nEvtRan)/((float)nEvtTotal);
  float trigevtfrac = ((float)nEvtTrigger)/((float)nEvtTotal);

  cout<<"---------------------------------------------"<<endl;
  cout<<"Summary:"<<endl;
  cout<<"Total events = "<<nEvtTotal<<endl;
  cout<<"Total events ran = "<<nEvtRan<<" ("<<goodevtfrac*100<<" %)"<<endl;
  cout<<"Total Triggered events = "<<nEvtTrigger<<" ("<<trigevtfrac*100<<" %)"<<endl;
  cout<<"---------------------------------------------"<<endl;
  cout<<"-----------------------------------------------"<<endl;
  cout<<"               Event Selection                 "<<endl;
  cout<<"-----------------------------------------------"<<endl;
  cout<<"Events with 4 or more Light Leptons = "<<n_l4<<endl;
  cout<<"Events with 3 Light Leptons         = "<<n_l3<<endl;
  cout<<"Events with 2LSS Light Leptons      = "<<n_l2SS<<endl;
  cout<<"Events with 1L2J Final state        = "<<n_l12j<<endl;
  
  time(&end);

  double time_taken = double(end-start);
  cout<<"Time taken by the program is= "<<fixed<<time_taken<<setprecision(5);
  cout<<"sec"<<endl;
}

void AnaScript::Terminate()
{
  
}

Bool_t AnaScript::Process(Long64_t entry)
{
  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC  .SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);

  //Setting verbosity:
  //Verbosity determines the number of processed events after which
  if(_verbosity==0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;

  nEvtTotal++;
  h.nevt->Fill(0);
  
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);

  GoodEvt = GoodEvt2016 && GoodEvt2017 && GoodEvt2018;

  if(GoodEvt){
    nEvtRan++; //only good events
    h.nevt->Fill(1);
    h.acceptance[0]->Fill(0);
    h.acceptance[1]->Fill(0);
    h.acceptance[2]->Fill(0);
    h.acceptance[3]->Fill(0);
    triggerRes=true; //Always true for MC
    
    if(_data==1){
      trigger2018 = (_year==2018 ? (_lep==1 ? *HLT_IsoMu24==1 : _lep==0 && *HLT_Ele32_WPTight_Gsf) : 1);
      trigger2017 = (_year==2017 ? (_lep==1 ? *HLT_IsoMu27==1 : _lep==0 && (*HLT_Ele32_WPTight_Gsf)) : 1);
      trigger2016 = (_year==2016 ? (_lep==1 ? (*HLT_IsoMu24==1) : _lep==0 && *HLT_Ele27_WPTight_Gsf) : 1);
           
      triggerRes = trigger2018 && trigger2017 && trigger2016;
     
    }
    if(triggerRes){
      nEvtTrigger++; //only triggered events
      h.nevt->Fill(2);

      /*******************************************************
       *               Gen particle block                    *
       *******************************************************/   
      //Clearing the Gen Arrays
      genMuon.clear();
      genElectron.clear();
      genLightLepton.clear();
      genVLLDoublet.clear();
      genVLLNuDoublet.clear();

      GenMuon();
      GenElectron();
      VLLDoublet();

      //Sorting the Gen Arrays
      SortPt(genMuon);
      SortPt(genElectron);
      SortPt(genVLLDoublet);
      SortPt(genVLLNuDoublet);
      
      //Make gen-level plots here.
      //VLL Doublet Plots      
      h.vllD[0]->Fill((int)genVLLDoublet.size());
      h.vllNuD[0]->Fill((int)genVLLNuDoublet.size());
      bool is_LL = false;
      bool is_LN = false;
      bool is_NN = false;
      if((int)genVLLDoublet.size()>1)
	if(abs(genVLLDoublet.at(0).pdgid) == 17 && abs(genVLLDoublet.at(1).pdgid) == 17){ //LL
	  h.decayMode[2]->Fill(1);
	  is_LL = true;
	}
      if((int)genVLLDoublet.size()>0 && (int)genVLLNuDoublet.size()>0)
	if(abs(genVLLDoublet.at(0).pdgid) == 17 && abs(genVLLNuDoublet.at(0).pdgid) == 18){ //LN
	  h.decayMode[2]->Fill(2);
	  is_LN = true;
	}
      if((int)genVLLNuDoublet.size()>1)
	if(abs(genVLLNuDoublet.at(0).pdgid) == 18 && abs(genVLLNuDoublet.at(1).pdgid) == 18){ //NN
	  h.decayMode[2]->Fill(3);
	  is_NN = true;
	}
      
      /*******************************************************
       *               Reco particle block                   *
       *******************************************************/   
      //Defining and Clearing the arrays
      Muon.clear();
      Electron.clear();
      LightLepton.clear();
      loosellep.clear();
      Tau.clear();
      Jets.clear();
      bJets.clear();

      //Calling the function that consists the electron and muons block
      createLightLeptons();
      createTaus();
      createJets();
      
      //Sorting the arrays
      SortPt(Muon);
      SortPt(Electron);
      SortPt(LightLepton);
      SortPt(loosellep);
      SortPt(Tau);
      SortPt(Jets);
      //Basic object-level plots:

      //Light Lepton Plots
      h.llep[0]->Fill((int)LightLepton.size());
      
      //ELectrons
      h.ele[0]->Fill((int)Electron.size());
      for(int i=0; i<(int)Electron.size(); i++){
	h.ele[1]->Fill(Electron.at(i).v.Pt());
	h.ele[2]->Fill(Electron.at(i).v.Eta());
	h.ele[3]->Fill(Electron.at(i).v.Phi());
      }
      if((int)Electron.size()>1){
	h.ele[4]->Fill(Electron.at(0).v.Pt());
	h.ele[5]->Fill(Electron.at(0).v.Eta());
	h.ele[6]->Fill(Electron.at(0).v.Phi());
	h.ele[7]->Fill(Electron.at(1).v.Pt());
	h.ele[8]->Fill(Electron.at(1).v.Eta());
	h.ele[9]->Fill(Electron.at(1).v.Phi());
	h.ele[10]->Fill(delta_phi(Electron.at(0).v.Phi(),*MET_phi));
	h.ele[11]->Fill(delta_phi(Electron.at(1).v.Phi(),*MET_phi));
      }
      
      //Muons
      h.mu[0]->Fill((int)Muon.size());
      for(int i=0; i<(int)Muon.size(); i++){
	h.mu[1]->Fill(Muon.at(i).v.Pt());
	h.mu[2]->Fill(Muon.at(i).v.Eta());
	h.mu[3]->Fill(Muon.at(i).v.Phi());
      }
      if((int)Muon.size()>1){
	h.mu[4]->Fill(Muon.at(0).v.Pt());
	h.mu[5]->Fill(Muon.at(0).v.Eta());
	h.mu[6]->Fill(Muon.at(0).v.Phi());
	h.mu[7]->Fill(Muon.at(1).v.Pt());
	h.mu[8]->Fill(Muon.at(1).v.Eta());
	h.mu[9]->Fill(Muon.at(1).v.Phi());
	h.mu[10]->Fill(delta_phi(Muon.at(0).v.Phi(),*MET_phi));
	h.mu[11]->Fill(delta_phi(Muon.at(1).v.Phi(),*MET_phi));
      }

      //MET
      h.met[0]->Fill(*MET_pt);
      h.met[1]->Fill(*MET_phi);

      //Filling Event level variable
      LT=0.0;
      LT_MET = 0.0;
      HT = 0.0;
      ST = 0.0;
      if((int)LightLepton.size()>0){
	for(int i=0; i<(int)LightLepton.size(); i++){
	  LT += LightLepton.at(i).v.Pt();
	}
	LT_MET = LT + *MET_pt;
	h.lt[0]->Fill(LT);
	h.lt_met[0]->Fill(LT_MET);
      }
      if((int)Jets.size()>1){
	for(int k=0; k<(int)Jets.size(); k++)
	  HT += Jets.at(k).v.Pt();
      }
      ST = LT_MET + HT;

      /***********************************************************************************
       *                                 Analysis block                                  *
       ***********************************************************************************/
      is_l2SS_event = false;
      is_l2OS_event = false;
      is_l12j_event = false;
      is_l3_event   = false;
      is_l4_event   = false;
      evt_trigger   = false;
      
      //4L Inclusive
      if((int)LightLepton.size() > 3){
	for(int i=0; i<4; i++){
	  if( fabs(LightLepton.at(i).id) == 11 && LightLepton.at(i).v.Pt() > 32 )
	    evt_trigger = true;
	  if( fabs(LightLepton.at(i).id) == 13 && LightLepton.at(i).v.Pt() > 26 )
	    evt_trigger = true;
	}
	if(evt_trigger){	  
	  is_l4_event = true;
	  n_l4++;
	  h.acceptance[0]->Fill(1);
	}
      }
      //Exactly 3L
      if((int)LightLepton.size() == 3){
	for(int i=0; i<3; i++){
	  if( fabs(LightLepton.at(i).id) == 11 && LightLepton.at(i).v.Pt() > 32 )
	    evt_trigger = true;
	  if( fabs(LightLepton.at(i).id) == 13 && LightLepton.at(i).v.Pt() > 26 )
	    evt_trigger = true;
	}
	if(evt_trigger){	  
	  is_l3_event = true;
	  n_l3++;
	  h.acceptance[0]->Fill(2);
	}
      }
      //2L SS
      if((int)LightLepton.size()==2){
	if((LightLepton.at(0).id)*(LightLepton.at(1).id) == 169 || (LightLepton.at(0).id)*(LightLepton.at(1).id) == 121 || (LightLepton.at(0).id)*(LightLepton.at(1).id) == 143){
	  for(int i=0; i<2; i++){
	    if( fabs(LightLepton.at(i).id) == 11 && LightLepton.at(i).v.Pt() > 32 )
	      evt_trigger = true;
	    if( fabs(LightLepton.at(i).id) == 13 && LightLepton.at(i).v.Pt() > 26 )
	      evt_trigger = true;
	  }
	  if(evt_trigger){	  
	    is_l2SS_event = true;
	    n_l2SS++;
	    h.acceptance[0]->Fill(3);
	  }
	}
      }
      //1L2J
      else if((int)LightLepton.size()==1 && (int)Jets.size()>1){
	if( fabs(LightLepton.at(0).id) == 11 && LightLepton.at(0).v.Pt() > 32 )
	  evt_trigger = true;
	if( fabs(LightLepton.at(0).id) == 13 && LightLepton.at(0).v.Pt() > 26 )
	  evt_trigger = true;
	if(evt_trigger){	  
	  is_l12j_event = true;
	  n_l12j++;
	  h.acceptance[0]->Fill(4);
	}
      }
      
      //Event level variables for 4L Events
      if(is_l4_event){
	//------>Kinematics
	h.hist_l4[0]->Fill((int)LightLepton.size());
	h.hist_l4[1]->Fill(LightLepton.at(0).v.Pt());
	h.hist_l4[2]->Fill(LightLepton.at(1).v.Pt());
	h.hist_l4[3]->Fill(LightLepton.at(2).v.Pt());
	h.hist_l4[4]->Fill(LightLepton.at(3).v.Pt());
	h.hist_l4[5]->Fill(LightLepton.at(0).v.Eta());
	h.hist_l4[6]->Fill(LightLepton.at(1).v.Eta());
	h.hist_l4[7]->Fill(LightLepton.at(2).v.Eta());
	h.hist_l4[8]->Fill(LightLepton.at(3).v.Eta());
	h.hist_l4[9]->Fill(LightLepton.at(0).v.Phi());
	h.hist_l4[10]->Fill(LightLepton.at(1).v.Phi());
	h.hist_l4[11]->Fill(LightLepton.at(2).v.Phi());
	h.hist_l4[12]->Fill(LightLepton.at(3).v.Phi());
	//------>Event Level Kinematics
	h.hist_l4[13]->Fill(LT);
	h.hist_l4[14]->Fill(*MET_pt);
	h.hist_l4[15]->Fill(LT_MET);
	//------>Angular Variables
	h.hist_l4[16]->Fill(delta_phi(LightLepton.at(0).v.Phi(),*MET_phi));
	h.hist_l4[17]->Fill(delta_phi(LightLepton.at(1).v.Phi(),*MET_phi));
	h.hist_l4[18]->Fill(delta_phi(LightLepton.at(2).v.Phi(),*MET_phi));
	h.hist_l4[19]->Fill(delta_phi(LightLepton.at(3).v.Phi(),*MET_phi));
	h.hist_l4[20]->Fill(delta_phi(LightLepton.at(0).v.Phi(),LightLepton.at(1).v.Phi()));
	h.hist_l4[21]->Fill(delta_phi(LightLepton.at(0).v.Phi(),LightLepton.at(2).v.Phi()));
	h.hist_l4[22]->Fill(delta_phi(LightLepton.at(0).v.Phi(),LightLepton.at(3).v.Phi()));
	h.hist_l4[23]->Fill(delta_phi(LightLepton.at(1).v.Phi(),LightLepton.at(2).v.Phi()));
	h.hist_l4[24]->Fill(delta_phi(LightLepton.at(1).v.Phi(),LightLepton.at(3).v.Phi()));
	h.hist_l4[25]->Fill(delta_phi(LightLepton.at(2).v.Phi(),LightLepton.at(3).v.Phi()));
	//------>Transverse Mass
	h.hist_l4[26]->Fill(transv_mass(LightLepton.at(0).v.Pt(),LightLepton.at(0).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l4[27]->Fill(transv_mass(LightLepton.at(1).v.Pt(),LightLepton.at(1).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l4[28]->Fill(transv_mass(LightLepton.at(2).v.Pt(),LightLepton.at(2).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l4[29]->Fill(transv_mass(LightLepton.at(3).v.Pt(),LightLepton.at(3).v.Phi(),*MET_pt,*MET_phi));
	//------>DeltaR Plots
	h.hist_l4[30]->Fill(LightLepton.at(0).v.DeltaR(LightLepton.at(1).v));
	h.hist_l4[31]->Fill(LightLepton.at(0).v.DeltaR(LightLepton.at(2).v));
	h.hist_l4[32]->Fill(LightLepton.at(0).v.DeltaR(LightLepton.at(3).v));
	h.hist_l4[33]->Fill(LightLepton.at(1).v.DeltaR(LightLepton.at(2).v));
	h.hist_l4[34]->Fill(LightLepton.at(1).v.DeltaR(LightLepton.at(3).v));
	h.hist_l4[35]->Fill(LightLepton.at(2).v.DeltaR(LightLepton.at(3).v));
	//Invariant Mass
	h.hist_l4[36]->Fill((LightLepton.at(0).v+LightLepton.at(1).v).M());
	h.hist_l4[37]->Fill((LightLepton.at(0).v+LightLepton.at(2).v).M());
	h.hist_l4[38]->Fill((LightLepton.at(0).v+LightLepton.at(3).v).M());
	h.hist_l4[39]->Fill((LightLepton.at(1).v+LightLepton.at(2).v).M());
	h.hist_l4[40]->Fill((LightLepton.at(1).v+LightLepton.at(3).v).M());
	h.hist_l4[41]->Fill((LightLepton.at(2).v+LightLepton.at(3).v).M());
      }
	
      //Event level Variables for 3L events
      if(is_l3_event){
	//------>Kinematics
	h.hist_l3[0]->Fill((int)LightLepton.size());
	h.hist_l3[1]->Fill(LightLepton.at(0).v.Pt());
	h.hist_l3[2]->Fill(LightLepton.at(1).v.Pt());
	h.hist_l3[3]->Fill(LightLepton.at(2).v.Pt());
	h.hist_l3[4]->Fill(LightLepton.at(0).v.Eta());
	h.hist_l3[5]->Fill(LightLepton.at(1).v.Eta());
	h.hist_l3[6]->Fill(LightLepton.at(2).v.Eta());
	h.hist_l3[7]->Fill(LightLepton.at(0).v.Phi());
	h.hist_l3[8]->Fill(LightLepton.at(1).v.Phi());
	h.hist_l3[9]->Fill(LightLepton.at(2).v.Phi());
	//------>Event Level Kinematics
	h.hist_l3[10]->Fill(LT);
	h.hist_l3[11]->Fill(*MET_pt);
	h.hist_l3[12]->Fill(LT_MET);
	//------>Angular Variables
	h.hist_l3[13]->Fill(delta_phi(LightLepton.at(0).v.Phi(),LightLepton.at(1).v.Phi()));
	h.hist_l3[14]->Fill(delta_phi(LightLepton.at(0).v.Phi(),LightLepton.at(2).v.Phi()));
	h.hist_l3[15]->Fill(delta_phi(LightLepton.at(1).v.Phi(),LightLepton.at(2).v.Phi()));
	h.hist_l3[16]->Fill(delta_phi(LightLepton.at(0).v.Phi(),*MET_phi));
	h.hist_l3[17]->Fill(delta_phi(LightLepton.at(1).v.Phi(),*MET_phi));
	h.hist_l3[18]->Fill(delta_phi(LightLepton.at(2).v.Phi(),*MET_phi));
	//------>Transverse Mass
	h.hist_l3[19]->Fill(transv_mass(LightLepton.at(0).v.Pt(),LightLepton.at(0).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l3[20]->Fill(transv_mass(LightLepton.at(1).v.Pt(),LightLepton.at(1).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l3[21]->Fill(transv_mass(LightLepton.at(2).v.Pt(),LightLepton.at(2).v.Phi(),*MET_pt,*MET_phi));
	//------>DeltaR Plots
	h.hist_l3[22]->Fill(LightLepton.at(0).v.DeltaR(LightLepton.at(1).v));
	h.hist_l3[23]->Fill(LightLepton.at(0).v.DeltaR(LightLepton.at(2).v));
	h.hist_l3[24]->Fill(LightLepton.at(1).v.DeltaR(LightLepton.at(2).v));
	//Invariant Mass
	h.hist_l3[25]->Fill((LightLepton.at(0).v+LightLepton.at(1).v).M());
	h.hist_l3[26]->Fill((LightLepton.at(0).v+LightLepton.at(2).v).M());
	h.hist_l3[27]->Fill((LightLepton.at(1).v+LightLepton.at(2).v).M());
      }

      //Event Level Variables for 2L SS events
      if(is_l2SS_event){
	//------>Kinematics
	h.hist_l2SS[0]->Fill((int)LightLepton.size());
	h.hist_l2SS[1]->Fill(LightLepton.at(0).v.Pt());
	h.hist_l2SS[2]->Fill(LightLepton.at(1).v.Pt());
	h.hist_l2SS[3]->Fill(LightLepton.at(0).v.Eta());
	h.hist_l2SS[4]->Fill(LightLepton.at(1).v.Eta());
	h.hist_l2SS[5]->Fill(LightLepton.at(0).v.Phi());
	h.hist_l2SS[6]->Fill(LightLepton.at(1).v.Phi());
	//------>Event Level Kinematics
	h.hist_l2SS[7]->Fill(LT);
	h.hist_l2SS[8]->Fill(*MET_pt);
	h.hist_l2SS[9]->Fill(LT_MET);
	//------>Angular Variables
	h.hist_l2SS[10]->Fill(delta_phi(LightLepton.at(0).v.Phi(),LightLepton.at(1).v.Phi()));
	h.hist_l2SS[11]->Fill(delta_phi(LightLepton.at(0).v.Phi(),*MET_phi));
	h.hist_l2SS[12]->Fill(delta_phi(LightLepton.at(1).v.Phi(),*MET_phi));
	//------>Transverse Mass
	h.hist_l2SS[13]->Fill(transv_mass(LightLepton.at(0).v.Pt(),LightLepton.at(0).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l2SS[14]->Fill(transv_mass(LightLepton.at(1).v.Pt(),LightLepton.at(1).v.Phi(),*MET_pt,*MET_phi));
	//------>DeltaR Plots
	h.hist_l2SS[15]->Fill(LightLepton.at(0).v.DeltaR(LightLepton.at(1).v));
	//Invariant Mass
	h.hist_l2SS[16]->Fill((LightLepton.at(0).v+LightLepton.at(1).v).M());
      }
      //Event level variables for 1L2J events
      if(is_l12j_event){
	//Kinematics
	h.hist_l12j[0]->Fill((int)LightLepton.size());
	h.hist_l12j[1]->Fill((int)Jets.size());
	h.hist_l12j[2]->Fill(LightLepton.at(0).v.Pt());
	h.hist_l12j[3]->Fill(Jets.at(0).v.Pt());
	h.hist_l12j[4]->Fill(Jets.at(1).v.Pt());
	h.hist_l12j[5]->Fill(LightLepton.at(0).v.Eta());
	h.hist_l12j[6]->Fill(Jets.at(0).v.Eta());
	h.hist_l12j[7]->Fill(Jets.at(1).v.Eta());
	h.hist_l12j[8]->Fill(LightLepton.at(0).v.Phi());
	h.hist_l12j[9]->Fill(Jets.at(0).v.Phi());
	h.hist_l12j[10]->Fill(Jets.at(1).v.Phi());
	//------>Transverse Mass
	h.hist_l12j[11]->Fill(transv_mass(LightLepton.at(0).v.Pt(),LightLepton.at(0).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l12j[12]->Fill(transv_mass(Jets.at(0).v.Pt(),Jets.at(0).v.Phi(),*MET_pt,*MET_phi));
	h.hist_l12j[13]->Fill(transv_mass(Jets.at(1).v.Pt(),Jets.at(1).v.Phi(),*MET_pt,*MET_phi));
	//Angular Variables
	h.hist_l12j[14]->Fill(delta_phi(LightLepton.at(0).v.Phi(),*MET_phi));
	h.hist_l12j[15]->Fill(delta_phi(Jets.at(0).v.Phi(),*MET_phi));
	h.hist_l12j[16]->Fill(delta_phi(Jets.at(1).v.Phi(),*MET_phi));
	h.hist_l12j[17]->Fill(delta_phi(Jets.at(0).v.Phi(),LightLepton.at(0).v.Phi()));
	h.hist_l12j[18]->Fill(delta_phi(Jets.at(1).v.Phi(),LightLepton.at(0).v.Phi()));
	h.hist_l12j[19]->Fill(delta_phi(Jets.at(0).v.Phi(),Jets.at(1).v.Phi()));
	h.hist_l12j[20]->Fill(delta_phi((Jets.at(0).v+Jets.at(1).v).Phi(),LightLepton.at(0).v.Phi()));
	h.hist_l12j[21]->Fill(delta_phi((Jets.at(0).v+Jets.at(1).v).Phi(),*MET_phi));
	h.hist_l12j[22]->Fill(delta_phi((Jets.at(0).v+Jets.at(1).v+LightLepton.at(0).v).Phi(),*MET_phi));
	h.hist_l12j[23]->Fill(delta_phi((Jets.at(0).v+Jets.at(1).v+LightLepton.at(0).v).Phi(),LightLepton.at(0).v.Phi()));
	h.hist_l12j[24]->Fill(delta_phi((Jets.at(0).v+Jets.at(1).v+LightLepton.at(0).v).Phi(),Jets.at(0).v.Phi()));
	h.hist_l12j[25]->Fill(delta_phi((Jets.at(0).v+Jets.at(1).v+LightLepton.at(0).v).Phi(),Jets.at(1).v.Phi()));
	//DiJet system
	h.hist_l12j[26]->Fill((Jets.at(0).v+Jets.at(1).v).M());
	h.hist_l12j[27]->Fill((Jets.at(0).v+Jets.at(1).v).Pt());
	h.hist_l12j[28]->Fill(transv_mass((Jets.at(0).v+Jets.at(1).v).Pt(),(Jets.at(0).v+Jets.at(1).v).Phi(),*MET_pt,*MET_phi));
	//dR
	h.hist_l12j[29]->Fill(Jets.at(0).v.DeltaR(Jets.at(1).v));
	//ljj system
	h.hist_l12j[30]->Fill((Jets.at(0).v+Jets.at(1).v+LightLepton.at(0).v).M());
	h.hist_l12j[31]->Fill((Jets.at(0).v+Jets.at(1).v+LightLepton.at(0).v).Pt());
	//Event Level Kinematics
	h.hist_l12j[32]->Fill(*MET_pt);
	h.hist_l12j[33]->Fill(HT);
	h.hist_l12j[34]->Fill(ST);
      }
      
    }//TriggeredEvts
  }//GoodEvt
  return kTRUE;
}
