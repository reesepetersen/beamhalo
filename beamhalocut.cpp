// Get some info out of a non-edm root file
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <numeric>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

using namespace std;

vector<string> split(string splitme, string splitter) {
  vector<size_t> splits;
  vector<string> splitted;
  size_t ifind = splitme.find(splitter);
  while (ifind != string::npos) {
    splits.push_back(ifind);
    ifind = splitme.find(splitter,ifind+1);
  }
  if (splits.size() > 0) {
    splitted.push_back(splitme.substr(0,splits[0]));
    for (int i=0; i<splits.size()-1; i++) {
      splitted.push_back(splitme.substr(splits[i]+splitter.length(),splits[i+1]-splits[i]-splitter.length()));
    }
    splitted.push_back(splitme.substr(splits[splits.size()-1]+splitter.length(),splitme.length()-splits[splits.size()-1]));
  }
  return splitted;
}

string barename(string filepath) {
  vector<string> pathvec = split(filepath,"/");
  string filename;
  if (pathvec.size() > 0) {
    filename = pathvec.back();
  }
  else {
    filename = filepath;
  }
  vector<string> fileextvec = split(filename,".");
  string bname = fileextvec[0];
  return bname;
}

int main(int argc, char** argv) {
    string infile = argv[1];
    TFile f(infile.c_str()); // Create TFile object for read file
    TTree* tree = (TTree*)f.Get("EventTree");
    tree->SetBranchStatus("*",1);
    string outname = barename(infile)+"_beamhalocut.root";
    TFile fnew(outname.c_str(),"RECREATE"); // Write over last output file of same name if it exists
    // new events tree
    TTree newtree("EventTree","Event data");
    Float_t pfMET;
    Float_t pfMETPhi;
    short nPho;
    vector<float>* phoE = nullptr;
    vector<float>* phoEt = nullptr;
    vector<float>* phoEta = nullptr;
    vector<float>* phoPhi = nullptr;
    vector<float>* phoR9 = nullptr;
    vector<float>* phoHoverE = nullptr;
    vector<float>* phoSigmaIEtaIEtaFull5x5 = nullptr;
    vector<float>* phoSigmaIEtaIPhiFull5x5 = nullptr;
    vector<float>* phoSigmaIPhiIPhiFull5x5 = nullptr;
    vector<float>* phoPFChWorstIso = nullptr;
    vector<float>* phoSCEta = nullptr;
    vector<float>* phoSCPhi = nullptr;
    vector<float>* phoSCE = nullptr;
    vector<float>* phoSCEtaWidth = nullptr;
    vector<float>* phoSCPhiWidth = nullptr;
    vector<float>* phoSeedTime = nullptr;
    Int_t nesRH;
    vector<float>* esrhE = nullptr;
    vector<float>* esrhiEta = nullptr;
    vector<float>* esrhiPhi = nullptr;
    vector<float>* esrhZside = nullptr;
    vector<float>* esrhPlane = nullptr;
    vector<float>* esrhStrip = nullptr;
    vector<float>* esrhEta = nullptr;
    vector<float>* esrhPhi = nullptr;
    vector<float>* esrhX = nullptr;
    vector<float>* esrhY = nullptr;
    vector<float>* esrhZ = nullptr;
    vector<float>* hbherhE = nullptr;
    vector<int>* hbherhDepth = nullptr;
    vector<int>* hbherhiEta = nullptr;
    vector<int>* hbherhiPhi = nullptr;
    vector<float>* hbherhEta = nullptr;
    vector<float>* hbherhPhi = nullptr;
    vector<float>* hbherhX = nullptr;
    vector<float>* hbherhY = nullptr;
    vector<float>* hbherhZ = nullptr;
    tree->SetBranchAddress("phoE",&phoE);
    tree->SetBranchAddress("phoEt",&phoEt);
    tree->SetBranchAddress("phoEta",&phoEta);
    tree->SetBranchAddress("phoPhi",&phoPhi);
    tree->SetBranchAddress("phoR9",&phoR9);
    tree->SetBranchAddress("phoHoverE",&phoHoverE);
    tree->SetBranchAddress("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5);
    tree->SetBranchAddress("phoSigmaIEtaIPhiFull5x5",&phoSigmaIEtaIPhiFull5x5);
    tree->SetBranchAddress("phoSigmaIPhiIPhiFull5x5",&phoSigmaIPhiIPhiFull5x5);
    tree->SetBranchAddress("phoPFChWorstIso",&phoPFChWorstIso);
    tree->SetBranchAddress("phoSCEta",&phoSCEta);
    tree->SetBranchAddress("phoSCPhi",&phoSCPhi);
    tree->SetBranchAddress("phoSCE",&phoSCE);
    tree->SetBranchAddress("phoSCEtaWidth",&phoSCEtaWidth);
    tree->SetBranchAddress("phoSCPhiWidth",&phoSCPhiWidth);
    tree->SetBranchAddress("phoSeedTime",&phoSeedTime);
    tree->SetBranchAddress("nesRH",&nesRH);
    tree->SetBranchAddress("esrhE",&esrhE);
    tree->SetBranchAddress("esrhiEta",&esrhiEta);
    tree->SetBranchAddress("esrhiPhi",&esrhiPhi);
    tree->SetBranchAddress("esrhZside",&esrhZside);
    tree->SetBranchAddress("esrhPlane",&esrhPlane);
    tree->SetBranchAddress("esrhStrip",&esrhStrip);
    tree->SetBranchAddress("esrhEta",&esrhEta);
    tree->SetBranchAddress("esrhPhi",&esrhPhi);
    tree->SetBranchAddress("esrhX",&esrhX);
    tree->SetBranchAddress("esrhY",&esrhY);
    tree->SetBranchAddress("esrhZ",&esrhZ);
    tree->SetBranchAddress("hbherhE",&hbherhE);
    tree->SetBranchAddress("hbherhDepth",&hbherhDepth);
    tree->SetBranchAddress("hbherhiEta",&hbherhiEta);
    tree->SetBranchAddress("hbherhiPhi",&hbherhiPhi);
    tree->SetBranchAddress("hbherhEta",&hbherhEta);
    tree->SetBranchAddress("hbherhPhi",&hbherhPhi);
    tree->SetBranchAddress("hbherhX",&hbherhX);
    tree->SetBranchAddress("hbherhY",&hbherhY);
    tree->SetBranchAddress("hbherhZ",&hbherhZ);
    float pho_total_count = 0.0f;
    int n = tree->GetEntries();
    float fn = n;
    vector<float> phoEnew;
    vector<float> phoEtnew;
    vector<float> phoEtanew;
    vector<float> phoPhinew;
    vector<float> phoR9new;
    vector<float> phoHoverEnew;
    vector<float> phoSigmaIEtaIEtaFull5x5new;
    vector<float> phoSigmaIEtaIPhiFull5x5new;
    vector<float> phoSigmaIPhiIPhiFull5x5new;
    vector<float> phoPFChWorstIsonew;
    vector<float> phoSCEtanew;
    vector<float> phoSCPhinew;
    vector<float> phoSCEnew;
    vector<float> phoSCEtaWidthnew;
    vector<float> phoSCPhiWidthnew;
    vector<float> phoSeedTimenew;
    Int_t nesRHnew;
    vector<float> esrhEnew;
    vector<float> esrhiEtanew;
    vector<float> esrhiPhinew;
    vector<float> esrhZsidenew;
    vector<float> esrhPlanenew;
    vector<float> esrhStripnew;
    vector<float> esrhEtanew;
    vector<float> esrhPhinew;
    vector<float> esrhXnew;
    vector<float> esrhYnew;
    vector<float> esrhZnew;
    vector<float> hbherhEnew;
    vector<int> hbherhDepthnew;
    vector<int> hbherhiEtanew;
    vector<int> hbherhiPhinew;
    vector<float> hbherhEtanew;
    vector<float> hbherhPhinew;
    vector<float> hbherhXnew;
    vector<float> hbherhYnew;
    vector<float> hbherhZnew;
    newtree.Branch("phoE",&phoEnew);
    newtree.Branch("phoEt",&phoEtnew);
    newtree.Branch("phoEta",&phoEtanew);
    newtree.Branch("phoPhi",&phoPhinew);
    newtree.Branch("phoR9",&phoR9new);
    newtree.Branch("phoHoverE",&phoHoverEnew);
    newtree.Branch("phoSigmaIEtaIEtaFull5x5",&phoSigmaIEtaIEtaFull5x5new);
    newtree.Branch("phoSigmaIEtaIPhiFull5x5",&phoSigmaIEtaIPhiFull5x5new);
    newtree.Branch("phoSigmaIPhiIPhiFull5x5",&phoSigmaIPhiIPhiFull5x5new);
    newtree.Branch("phoPFChWorstIso",&phoPFChWorstIsonew);
    newtree.Branch("phoSCEta",&phoSCEtanew);
    newtree.Branch("phoSCPhi",&phoSCPhinew);
    newtree.Branch("phoSCE",&phoSCEnew);
    newtree.Branch("phoSCEtaWidth",&phoSCEtaWidthnew);
    newtree.Branch("phoSCPhiWidth",&phoSCPhiWidthnew);
    newtree.Branch("phoSeedTime",&phoSeedTimenew);
    newtree.Branch("nesRH",&nesRHnew);
    newtree.Branch("esrhE",&esrhEnew);
    newtree.Branch("esrhiEta",&esrhiEtanew);
    newtree.Branch("esrhiPhi",&esrhiPhinew);
    newtree.Branch("esrhZside",&esrhZsidenew);
    newtree.Branch("esrhPlane",&esrhPlanenew);
    newtree.Branch("esrhStrip",&esrhStripnew);
    newtree.Branch("esrhEta",&esrhEtanew);
    newtree.Branch("esrhPhi",&esrhPhinew);
    newtree.Branch("esrhX",&esrhXnew);
    newtree.Branch("esrhY",&esrhYnew);
    newtree.Branch("esrhZ",&esrhZnew);
    newtree.Branch("hbherhE",&hbherhEnew);
    newtree.Branch("hbherhDepth",&hbherhDepthnew);
    newtree.Branch("hbherhiEta",&hbherhiEtanew);
    newtree.Branch("hbherhiPhi",&hbherhiPhinew);
    newtree.Branch("hbherhEta",&hbherhEtanew);
    newtree.Branch("hbherhPhi",&hbherhPhinew);
    newtree.Branch("hbherhX",&hbherhXnew);
    newtree.Branch("hbherhY",&hbherhYnew);
    newtree.Branch("hbherhZ",&hbherhZnew);
    float phoEta_pass_count = 0.0f;
    float phoPhi_pass_count = 0.0f;
    float phoEt_pass_count = 0.0f;
    for (Int_t i=0; i<n; i++) {
      tree->GetEntry(i);
      for (int j=0; j<phoEta->size(); j++) {
        pho_total_count++;
        if(fabs(phoEta->at(j)) > 1.566 && fabs(phoEta->at(j)) < 2.5) {
          phoEta_pass_count++;
          if (fabs(phoPhi->at(j)) <= 0.2 || fabs(phoPhi->at(j)) >= 2.9416) {
            phoPhi_pass_count++;
            if (phoEt->at(j) > 250) {
              phoEt_pass_count++;
              phoEnew.push_back(phoE->at(j));
              phoEtnew.push_back(phoEt->at(j));
              phoEtanew.push_back(phoEta->at(j));
              phoPhinew.push_back(phoPhi->at(j));
              phoR9new.push_back(phoR9->at(j));
              phoHoverEnew.push_back(phoHoverE->at(j));
              phoSigmaIEtaIEtaFull5x5new.push_back(phoSigmaIEtaIEtaFull5x5->at(j));
              phoSigmaIEtaIPhiFull5x5new.push_back(phoSigmaIEtaIPhiFull5x5->at(j));
              phoSigmaIPhiIPhiFull5x5new.push_back(phoSigmaIPhiIPhiFull5x5->at(j));
              phoPFChWorstIsonew.push_back(phoPFChWorstIso->at(j));
              phoSCEtanew.push_back(phoSCEta->at(j));
              phoSCPhinew.push_back(phoSCPhi->at(j));
              phoSCEnew.push_back(phoSCE->at(j));
              phoSCEtaWidthnew.push_back(phoSCEtaWidth->at(j));
              phoSCPhiWidthnew.push_back(phoSCPhiWidth->at(j));
              phoSeedTimenew.push_back(phoSeedTime->at(j));
            }
          }
        }
      }
      if (phoEnew.size() >= 1) {
        nesRHnew = nesRH;
        esrhEnew = *esrhE;
        esrhiEtanew = *esrhiEta;
        esrhiPhinew = *esrhiPhi;
        esrhZsidenew = *esrhZside;
        esrhPlanenew = *esrhPlane;
        esrhStripnew = *esrhStrip;
        esrhEtanew = *esrhEta;
        esrhPhinew = *esrhPhi;
        esrhXnew = *esrhX;
        esrhYnew = *esrhY;
        esrhZnew = *esrhZ;
        hbherhEnew = *hbherhE;
        hbherhDepthnew = *hbherhDepth;
        hbherhiEtanew = *hbherhiEta;
        hbherhiPhinew = *hbherhiPhi;
        hbherhEtanew = *hbherhEta;
        hbherhPhinew = *hbherhPhi;
        hbherhXnew = *hbherhX;
        hbherhYnew = *hbherhY;
        hbherhZnew = *hbherhZ;
        newtree.Fill();
      }
      phoEnew.clear();
      phoEtnew.clear();
      phoEtanew.clear();
      phoPhinew.clear();
      phoR9new.clear();
      phoHoverEnew.clear();
      phoSigmaIEtaIEtaFull5x5new.clear();
      phoSigmaIEtaIPhiFull5x5new.clear();
      phoSigmaIPhiIPhiFull5x5new.clear();
      phoPFChWorstIsonew.clear();
      phoSCEtanew.clear();
      phoSCPhinew.clear();
      phoSCEnew.clear();
      phoSCEtaWidthnew.clear();
      phoSCPhiWidthnew.clear();
      phoSeedTimenew.clear();
      esrhEnew.clear();
      esrhiEtanew.clear();
      esrhiPhinew.clear();
      esrhZsidenew.clear();
      esrhPlanenew.clear();
      esrhStripnew.clear();
      esrhEtanew.clear();
      esrhPhinew.clear();
      esrhXnew.clear();
      esrhYnew.clear();
      esrhZnew.clear();
      hbherhEnew.clear();
      hbherhDepthnew.clear();
      hbherhiEtanew.clear();
      hbherhiPhinew.clear();
      hbherhEtanew.clear();
      hbherhPhinew.clear();
      hbherhXnew.clear();
      hbherhYnew.clear();
      hbherhZnew.clear();
    }
    newtree.Write();
    // new stats tree
    TTree* statstree = (TTree*)f.Get("statstree");
    statstree->SetBranchStatus("*",1);
    TTree newstatstree("statstree","Stats_Tree");
    Float_t n_events_precut;
    Float_t n_events_t_postcut;
    Float_t n_events_tm_postcut;
    statstree->SetBranchAddress("n_events_precut",&n_events_precut);
    statstree->SetBranchAddress("n_events_t_postcut",&n_events_t_postcut);
    statstree->SetBranchAddress("n_events_tm_postcut",&n_events_tm_postcut);
    Float_t n_events_precutnew;
    Float_t n_events_t_postcutnew;
    Float_t n_events_tm_postcutnew;
    Float_t n_events_eta_postcut;
    Float_t n_events_et_postcut;
    Float_t n_events_phi_postcut;
    newstatstree.Branch("n_events_precut",&n_events_precutnew,"n_events_precut/F");
    newstatstree.Branch("n_events_t_postcut",&n_events_t_postcutnew,"n_events_t_postcut/F");
    newstatstree.Branch("n_events_tm_postcut",&n_events_tm_postcutnew,"n_events_tm_postcut/F");
    newstatstree.Branch("n_events_eta_postcut",&n_events_eta_postcut,"n_events_eta_postcut/F");
    newstatstree.Branch("n_events_et_postcut",&n_events_et_postcut,"n_events_et_postcut/F");
    newstatstree.Branch("n_events_phi_postcut",&n_events_phi_postcut,"n_events_phi_postcut/F");
    n_events_precutnew = n_events_precut;
    n_events_t_postcutnew = n_events_t_postcut;
    n_events_tm_postcutnew = n_events_tm_postcut;
    n_events_eta_postcut = phoEta_pass_count;
    n_events_et_postcut = phoEt_pass_count;
    n_events_phi_postcut = phoPhi_pass_count;
    newstatstree.Fill();
    newstatstree.Write();
    fnew.Close();
    f.Close();
  return 0;
}
