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
  string infile = argv[1]; // Save input as string
  TFile f (infile.c_str()); // Create TFile object for read file
  TTree* tree = (TTree*)f.Get("EventTree");
  tree->SetBranchStatus("*",1);
  Long64_t event;
  vector<float>* phoE = nullptr;
  vector<float>* phoEt = nullptr;
  vector<float>* phoEta = nullptr;
  vector<float>* phoPhi = nullptr;
  vector<float>* phoSCEta = nullptr;
  vector<float>* phoSCPhi = nullptr;
  vector<float>* phoSCE = nullptr;
  vector<float>* phoSCPhiWidth = nullptr;
  vector<float>* phoSCEtaWidth = nullptr;
  vector<float>* eleEta = nullptr;
  vector<float>* elePhi = nullptr;
  vector<float>* AK4CHSJet_Eta = nullptr;
  vector<float>* AK4CHSJet_Phi = nullptr;
  vector<float>* muEta = nullptr;
  vector<float>* muPhi = nullptr;
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
  vector<int>* esrhHaloPreID = nullptr;
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("phoE",&phoE);
  tree->SetBranchAddress("phoEt",&phoEt);
  tree->SetBranchAddress("phoEta",&phoEta);
  tree->SetBranchAddress("phoPhi",&phoPhi);
  tree->SetBranchAddress("phoSCEta",&phoSCEta);
  tree->SetBranchAddress("phoSCPhi",&phoSCPhi);
  tree->SetBranchAddress("phoSCE",&phoSCE);
  tree->SetBranchAddress("phoSCEtaWidth",&phoSCEtaWidth);
  tree->SetBranchAddress("phoSCPhiWidth",&phoSCPhiWidth);
  tree->SetBranchAddress("eleEta",&eleEta);
  tree->SetBranchAddress("elePhi",&elePhi);
  tree->SetBranchAddress("AK4CHSJet_Eta",&AK4CHSJet_Eta);
  tree->SetBranchAddress("AK4CHSJet_Phi",&AK4CHSJet_Phi);
  tree->SetBranchAddress("muEta",&muEta);
  tree->SetBranchAddress("muPhi",&muPhi);
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
  tree->SetBranchAddress("esrhHaloPreID",&esrhHaloPreID);
  int nesrh;
  float ecal_z = 319.5;//cm
  float obR;
  float obphi;
  float obeta;
  float obx;
  float oby;
  cout << "esrhX,esrhY\n";
  //int n = tree->GetEntriesFast();
  //for (Int_t i=0; i<n; i++) {
    //tree->GetEntry(i);
    tree->GetEntry(0);
    nesrh = esrhE->size();
    for(int j = 0; j < nesrh; j++) {
      cout << esrhX->at(j) << "," << esrhY->at(j) << '\n';
    }
  cout << "phoX,phoY\n";
  int npho = phoEta->size();
  for (int i = 0; i < npho; i++) {
    phophi = phoPhi->at(i);
    phoeta = phoEta->at(i);
    phoR = ecal_z/sinh(phoeta);
    phox = phoR*cos(phophi);
    phoy = phoR*sin(phophi);
    cout << phox << "," << phoy << '\n';
  }
  //}
  f.Close();
  return 0;
}
