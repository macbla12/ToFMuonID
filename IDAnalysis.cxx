
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <string>
#include <TLegend.h>
#include <vector>
#include <tuple>

#include "ToFSim.cxx"



void IDAnalysis()
{
   //////////////////////
   //Setting up constants
   //////////////////////

   static double MuonMass=0.1056583;
   static double ElectronMass=0.00051099895;
   static double PionMass=0.13957039;
   static double KaonMass=0.493677;
   static double ProtonMass=0.93827208;

   static double c_light=299.792458;

   vector<double> masses = {MuonMass, ElectronMass, PionMass, KaonMass, ProtonMass};
   vector<int> pdgs = {13, 11, 211, 321, 2212};

   gROOT->SetBatch(kTRUE);
   gROOT->ProcessLine("gErrorIgnoreLevel = 3000;");
   //gStyle->SetOptStat(0);

   double DEG=180/TMath::Pi();

   TFile *PaperFile= new TFile("Plots/ToFPlots.root", "RECREATE");


   //////////////////////
   //Setting up histograms
   //////////////////////
   static constexpr int NumOfFiles=5;
   TH1D *BToFEtaRangeHist[NumOfFiles], *BToFPhiRangeHist[NumOfFiles], *EToFEtaRangeHist[NumOfFiles], *EToFPhiRangeHist[NumOfFiles];
   TH1D *ToFBTimeHist[NumOfFiles], *ToFETimeHist[NumOfFiles];
   TH1D *NumberOfBHits[NumOfFiles], *NumberOfEHits[NumOfFiles];
   TH1D *ToFTimeHist[NumOfFiles], *ToFMassHist[NumOfFiles];
   TH1D *DoubleToFMassHist[NumOfFiles];
   TH1D *DistanceHist= new TH1D("DistanceHist","DistanceHist",50,0,1000);
   TH1D *LengthDiffrance= new TH1D("LengthDiffrance","LengthDiffrance",50,-100,2000);
   TH2D *BetaVsMom= new TH2D("BetaVsMom","BetaVsMom",100,0,5,100,0,1.1);
   
   vector<TString> files(NumOfFiles);
   
   //files.at(0)="/run/media/epic/Data/Muons/Grape-10x275/Current/reco*.root";
   //files.at(0)="/run/media/epic/Data/Tau/reco/Energy_10x275/double_pi/recoDoublePi.root";

   files.at(0)="/run/media/epic/Data/Background/Muons/Continuous/reco_*.root";
   files.at(1)="/run/media/epic/Data/Background/Pions/Continuous/reco_*.root";
   files.at(2)="/run/media/epic/Data/Background/SingleParticles/SingleFiles/Electrons.root";
   files.at(3)="/run/media/epic/Data/Background/SingleParticles/SingleFiles/Kaons.root";
   files.at(4)="/run/media/epic/Data/Background/SingleParticles/SingleFiles/Protons.root";

   
   for(int File=0; File<NumOfFiles;File++)
   {
      string name;
      if(File==0) name="Muons";
      if(File==1) name="Pions";
      if(File==2) name="Electrons";
      if(File==3) name="Kaons";
      if(File==4) name="Protons";

      // Set up input file chain
      TChain *mychain = new TChain("events");
      mychain->Add(files.at(File));

      // Initialize reader
      TTreeReader tree_reader(mychain);

      // Get Particle Information
      TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
      TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
      TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
      TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
      TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
      TTreeReaderArray<double> partMass(tree_reader, "MCParticles.mass");
      TTreeReaderArray<float> partCharge(tree_reader, "MCParticles.charge");
      TTreeReaderArray<unsigned int> partParb(tree_reader, "MCParticles.parents_begin");
      TTreeReaderArray<unsigned int> partPare(tree_reader, "MCParticles.parents_end");
      TTreeReaderArray<int> partParI(tree_reader, "_MCParticles_parents.index");

      // Get Reconstructed Track Information
      TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
      TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
      TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
      TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
      TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
      TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
      TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");
      TTreeReaderArray<float> trackRefX(tree_reader, "ReconstructedChargedParticles.referencePoint.x");
      TTreeReaderArray<float> trackRefY(tree_reader, "ReconstructedChargedParticles.referencePoint.y");
      TTreeReaderArray<float> trackRefZ(tree_reader, "ReconstructedChargedParticles.referencePoint.z");


      // Get Reconstructed ToF Barrel Information
      TTreeReaderArray<float> BToFPosX(tree_reader, "TOFBarrelRecHits.position.x");
      TTreeReaderArray<float> BToFPosY(tree_reader, "TOFBarrelRecHits.position.y");
      TTreeReaderArray<float> BToFPosZ(tree_reader, "TOFBarrelRecHits.position.z");
      TTreeReaderArray<float> BToFTime(tree_reader, "TOFBarrelRecHits.time");

      // Get Reconstructed ToF Endcap Information
      TTreeReaderArray<float> EToFPosX(tree_reader, "TOFEndcapRecHits.position.x");
      TTreeReaderArray<float> EToFPosY(tree_reader, "TOFEndcapRecHits.position.y");
      TTreeReaderArray<float> EToFPosZ(tree_reader, "TOFEndcapRecHits.position.z");
      TTreeReaderArray<float> EToFTime(tree_reader, "TOFEndcapRecHits.time");


     
      //==================================//
      BToFEtaRangeHist[File]= new TH1D(Form("BToFEtaRangeHist%s",name.c_str()),Form("BToFEtaRangeHist%s",name.c_str()),50,-2.5,2.5);
      BToFPhiRangeHist[File]= new TH1D(Form("BToFPhiRangeHist%s",name.c_str()),Form("BToFPhiRangeHist%s",name.c_str()),50,-180,180);

      EToFEtaRangeHist[File]= new TH1D(Form("EToFEtaRangeHist%s",name.c_str()),Form("EToFEtaRangeHist%s",name.c_str()),50,1.5,3.5);
      EToFPhiRangeHist[File]= new TH1D(Form("EToFPhiRangeHist%s",name.c_str()),Form("EToFPhiRangeHist%s",name.c_str()),50,-180,180);

      ToFBTimeHist[File]= new TH1D(Form("ToFBarrelTimeHist%s",name.c_str()),Form("ToFBarrelTimeHist%s",name.c_str()),50,0,15);
      ToFETimeHist[File]= new TH1D(Form("ToFEndcapTimeHist%s",name.c_str()),Form("ToFEndcapTimeHist%s",name.c_str()),50,0,15);

      NumberOfBHits[File]= new TH1D(Form("NumberOfBarrelHits%s",name.c_str()),Form("NumberOfBarrelHits%s",name.c_str()),10,-0.5,9.5);
      NumberOfEHits[File]= new TH1D(Form("NumberOfEndcapHits%s",name.c_str()),Form("NumberOfEndcapHits%s",name.c_str()),10,-0.5,9.5);
      




      ToFMassHist[File]= new TH1D(Form("ToFMassHist%s",name.c_str()),Form("ToFMassHist%s",name.c_str()),50,0,0.3);
      //==================================//
      DoubleToFMassHist[File]= new TH1D(Form("DoubleToFMassHist%s",name.c_str()),Form("ToFMassHist%s",name.c_str()),50,0,0.3);

      //====================================//
     



      int eventID=0;
      double FoundParticles=0;
      double particscount=0;
      double goodcount=0;
      double badcount=0;
      double sigma_t_Barrel= 0.025;
      double sigma_t_Endcap= 0.04;


      while(tree_reader.Next()){
         eventID++;
         //if(eventID>40) break;
         if(eventID>20000) break;
         
         //if(File==0) if(eventID>3800) break;
         if(eventID%20000==0) cout<<"File "<<name<<" and event number... "<<eventID<<endl;

         //if(File==0) continue;
         //if(File==1) continue;
         //if(File==2) continue;
         //if(File==3) continue;
         //if(File==4) continue;
         double FilePDG;
         if(File==0) FilePDG=13;
         if(File==1) FilePDG=211;
         if(File==2) FilePDG=11;
         if(File==3) FilePDG=321;
         if(File==4) FilePDG=2212;

         
         int NBarrelHits=0;
         int NEndcapHits=0;
         for(int particle=0; particle<trackEng.GetSize(); particle++)
         {
            //Obligatory Cuts 
            double mass;

            int Found=0;
            TLorentzVector Partic;
            Partic.SetPxPyPzE(trackMomX[particle],trackMomY[particle],trackMomZ[particle],trackEng[particle]);
            if(Partic.P()>5) continue;
            particscount++;
            
            for(int hits=0; hits<BToFTime.GetSize(); hits++)
            {
               TVector3 ToFPos(BToFPosX[hits],BToFPosY[hits],BToFPosZ[hits]);
               Found=1;

               BToFEtaRangeHist[File]->Fill(ToFPos.Eta());
               BToFPhiRangeHist[File]->Fill(ToFPos.Phi()*DEG);
               
               ToFBTimeHist[File]->Fill(BToFTime[hits]);
               NBarrelHits++;
               auto [length, DistanceCheck, smallest_distance] = ToFSim(Partic, trackCharge[particle], ToFPos);
               DistanceHist->Fill(smallest_distance);
               LengthDiffrance->Fill(ToFPos.Mag() - length);

               //cout<<"Distance to BToF: "<<smallest_distance<<" cm, magnitude: "<<ToFPos.Mag()<<" length "<<length<<" cm"<<endl;
               if(!DistanceCheck) continue;
               
               double Beta = length/BToFTime[hits]/c_light;

               BetaVsMom->Fill(Partic.P(),Beta);
               int bestPDG = 0;
               double bestNSigma = 1e9;
               double nSigma[5];  
               for (int i = 0; i < pdgs.size(); i++) {

                  double t_exp = computeExpectedTime(masses[i], Partic.P(), length);
                  double dt = BToFTime[hits] - t_exp;

                  double n = fabs(dt) / sigma_t_Barrel;
                  nSigma[i] = n;

                  if (n < bestNSigma) {
                     bestNSigma = n;
                     bestPDG = pdgs[i];
                  }
               }


               if (bestNSigma > 3.0)
                  continue; 

            
               if(bestPDG==FilePDG) goodcount++;
               else badcount++;
               

               
            }
            for(int hits=0; hits<EToFTime.GetSize(); hits++)
            {
               TVector3 ToFPos(EToFPosX[hits],EToFPosY[hits],EToFPosZ[hits]);
               Found=1;

               EToFEtaRangeHist[File]->Fill(ToFPos.Eta());
               EToFPhiRangeHist[File]->Fill(ToFPos.Phi()*DEG);
               
               ToFETimeHist[File]->Fill(EToFTime[hits]);
               NEndcapHits++;
               auto [length, DistanceCheck, smallest_distance] = ToFSim(Partic, trackCharge[particle], ToFPos);
               DistanceHist->Fill(smallest_distance);
               LengthDiffrance->Fill(ToFPos.Mag() - length);
               if(!DistanceCheck) continue;

               //cout<<"Distance to EToF: "<<smallest_distance<<" cm, magnitude: "<<ToFPos.Mag()<<" length "<<length<<" cm"<<endl;
               double Beta = length/EToFTime[hits]/c_light;
               BetaVsMom->Fill(Partic.P(),Beta);
               int bestPDG = 0;
               double bestNSigma = 1e9;
               double nSigma[5];  
               for (int i = 0; i < pdgs.size(); i++) {

                  double t_exp = computeExpectedTime(masses[i], Partic.P(), length);
                  double dt = EToFTime[hits] - t_exp;

                  double n = fabs(dt) / sigma_t_Endcap;
                  nSigma[i] = n;

                  if (n < bestNSigma) {
                     bestNSigma = n;
                     bestPDG = pdgs[i];
                  }
               }


               if (bestNSigma > 3.0)
                  continue;  

 
               if(FilePDG==bestPDG) goodcount++;
               else badcount++;
               


            } 
            FoundParticles+=Found;
            
         }
         NumberOfBHits[File]->Fill(NBarrelHits);
         NumberOfEHits[File]->Fill(NEndcapHits);
         
      }
      

      cout<<"==========================="<<endl;
      cout<<"End of "<< name << " file"<<endl;
      cout<<"Number of events: "<<eventID<<endl;
      cout<<"Found particles: "<<FoundParticles<<"   All particles: "<<particscount<<endl;
      cout<<"Found particles: "<<FoundParticles<<endl;

      cout<<"Found Ratio: "<<FoundParticles*100/particscount<<'%'<<endl;
      cout<<"Purity: "<<goodcount*100/(goodcount+badcount)<<'%'<<endl;

      cout<<"==========================="<<endl;


      BToFEtaRangeHist[File]->Write();
      BToFPhiRangeHist[File]->Write();


      EToFEtaRangeHist[File]->Write();
      EToFPhiRangeHist[File]->Write();

      ToFBTimeHist[File]->Write();
      ToFETimeHist[File]->Write();

      NumberOfBHits[File]->Write();
      NumberOfEHits[File]->Write();
      
      ToFMassHist[File]->Write();
      DoubleToFMassHist[File]->Write();
   }
   
   DistanceHist->Write();
   LengthDiffrance->Write();
   BetaVsMom->Write();
   PaperFile->Close();


 
}

