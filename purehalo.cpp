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
    string fileout = barename(infile)+"_ph.root";
    TFile fnew (fileout.c_str(),"RECREATE"); // Write over last output file of same name if it exists
    tree->SetBranchStatus("*",1);
    auto newtree = tree->CloneTree(0);
    ULong64_t HLTPho;
    Float_t pfMET;
    UShort_t beamHaloSummary;
    vector<float>* phoEta = nullptr;
    vector<float>* phoE = nullptr;
    vector<float>* phoEt = nullptr; 
    tree->SetBranchAddress("HLTPho",&HLTPho);
    tree->SetBranchAddress("pfMET",&pfMET);
    tree->SetBranchAddress("beamHaloSummary",&beamHaloSummary);
    tree->SetBranchAddress("phoEta",&phoEta);
    tree->SetBranchAddress("phoEt",&phoEt);
    tree->SetBranchAddress("phoE",&phoE);
    int n = tree->GetEntriesFast();
    int npho;
    float phoeta;
    bool phopass;
    cout << "got to for loop\n";
    for (int i=0; i<n; i++) {
      cout << "start event loop\n";
      tree->GetEntry(i);
      cout << "got entry\n";
      if (HLTPho>>9&1 == 1 and pfMET > 200.0f and beamHaloSummary>>2&1==1) {
        cout << "pass event level\n";
        phopass = false;
        npho = phoEta->size();
        for (int j = 0; j < npho; j++) {
          cout << "start photon loop\n";
          phoeta = phoEta->at(j);
          if (fabs(phoeta) >= 1.65 and fabs(phoeta) <= 1.8 and phoE->at(j) < 5000 and phoEt->at(j) > 50) phopass = true;
        }
        if (phopass) newtree->Fill();
      }
    }
    newtree->Write();
    fnew.Close(); // Close files
    f.Close();
  return 0;
}
