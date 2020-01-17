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
    TFile *f = new TFile(infile.c_str()); // Create TFile object for read file
    auto *ggN = f->Get("ggNtuplizer"); // Get the TDirectoryFile object in the file called "ggNtuplizer"
    auto *ggN2 = dynamic_cast<TDirectoryFile*> (ggN); // It came back as a TObject, so now lets recast it as a TDirectoryFile so we can get the TTree inside
    auto *etree = ggN2->Get("EventTree"); // Get the TTree called EventTree, if there are multiple TTrees get the one with the lowest cycle number
    string outpath = "";
    string outname = barename(infile)+"_temcut.root";
    string fileout = outpath+outname;
    TFile *fnew = new TFile(fileout.c_str(),"RECREATE"); // Write over last output file of same name if it exists
    auto *tree = dynamic_cast<TTree*> (etree);
    tree->SetBranchStatus("*",1);
    auto newtree = tree->CloneTree(0);
    TTree* statstree = new TTree("statstree","Stats_Tree");
    ULong64_t HLTPho; // HLTPho
    Float_t pfMET; // pfMET
    tree->SetBranchAddress("HLTPho",&HLTPho);
    tree->SetBranchAddress("pfMET",&pfMET);
    float HLTPho_pass_count = 0.0f;
    float pfMET_pass_count = 0.0f;
    int n = tree->GetEntriesFast();
    float fn = n;
    for (Int_t i=0; i<n; i++) {
      tree->GetEntry(i);
      if (HLTPho>>9&1 == 1) {
        HLTPho_pass_count+=1.0f;
        if (pfMET > 150.0f) {
          pfMET_pass_count+=1.0f;
            newtree->Fill();
        }
      }
    }
    Float_t n_events_precut;
    statstree->Branch("n_events_precut",&n_events_precut,"n_events_precut/F");
    n_events_precut = fn;
    Float_t n_events_t_postcut;
    statstree->Branch("n_events_t_postcut",&n_events_t_postcut,"n_events_t_postcut/F");
    n_events_t_postcut = HLTPho_pass_count;
    Float_t n_events_tm_postcut;
    statstree->Branch("n_events_tm_postcut",&n_events_tm_postcut,"n_events_tm_postcut/F");
    n_events_tm_postcut = pfMET_pass_count;
    statstree->Fill();
    newtree->Write();
    statstree->Write();
    fnew->Close(); // Close files
    f->Close();
  return 0;
}
