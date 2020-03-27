// Get some info out of a non-edm root file
#include <cstdio>
#include <cstdlib>
#include <random>
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

bool check_samesign_eta(float phoEta, float hbherhEta) {
  if (phoEta > 0 && hbherhEta > 0) return true;
  else if (phoEta < 0 && hbherhEta < 0) return true;
  else return false;
}

int get_close_phi_herh(float phoPhi,vector<float> *hbherhPhi,vector<float> *hbherhEta) {
  int herh_count = 0;
  float delphi;
  float pi = 3.1415926535897932;
  for (int i = 0; i < hbherhPhi->size(); i++) {
    if (fabs(hbherhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
    delphi = fabs(hbherhPhi->at(i)-phoPhi);
    if (delphi > pi) delphi = 2*pi - delphi;
    if (delphi < 0.09) herh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
  }
  return herh_count;
}
            
int get_close_phi_herhzside(float phoPhi, float phoEta, vector<float> *hbherhPhi, vector<float> *hbherhEta) {
  int herh_count = 0;
  float delphi;
  float pi = 3.1415926535897932;
  for (int i = 0; i < hbherhPhi->size(); i++) {
    if (fabs(hbherhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
    delphi = fabs(hbherhPhi->at(i)-phoPhi);
    if (delphi > pi) delphi = 2*pi - delphi;
    bool samesign_eta = check_samesign_eta(phoEta,hbherhEta->at(i));
    if (delphi < 0.09 && samesign_eta) herh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
  }
  return herh_count;
}
            
int get_close_phi_esrh(float phoPhi,vector<float> *esrhPhi,vector<float> *esrhEta) {
  int esrh_count = 0;
  float delphi;
  float pi = 3.1415926535897932;
  for (int i = 0; i < esrhPhi->size(); i++) {
    if (fabs(esrhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
    delphi = fabs(esrhPhi->at(i)-phoPhi);
    if (delphi > pi) delphi = 2*pi - delphi;
    if (delphi < 0.09) esrh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
  }
  return esrh_count;
}

int get_close_phi_esrhzside(float phoPhi, float phoEta, vector<float> *esrhPhi,vector<float> *esrhEta) {
  int esrh_count = 0;
  float delphi;
  float pi = 3.1415926535897932;
  for (int i = 0; i < esrhPhi->size(); i++) {
    if (fabs(esrhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
    delphi = fabs(esrhPhi->at(i)-phoPhi);
    if (delphi > pi) delphi = 2*pi - delphi;
    bool samesign_eta = check_samesign_eta(phoEta,esrhEta->at(i));
    if (delphi < 0.09 && samesign_eta) esrh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
  }
  return esrh_count;
}

vector<float> phoAHETotal(vector<float> *esrhX, vector<float> *esrhY, vector<float> *esrhZ, vector<float> *esrhPhi, vector<float> *esrhE, vector<float> *esrhEta, vector<float> *hbherhX, vector<float> *hbherhY, vector<float> *hbherhZ, vector<float> *hbherhPhi, vector<float> *hbherhE, vector<float> *hbherhEta, vector<float> *phoEta, vector<float> *phoPhi) {
  vector<float> TotalVec;
  float pi = 3.14159265358979323846;
  float ecal_z = 319.5;//cm
  float alpha = pi/180;//1 degree in radians
  float phoR;
  float phophi;
  float esrhR;
  float esrhphi;
  float hbherhR;
  float hbherhphi;
  float delphi;
  float delr;
  float esrhEtot;
  float hbherhEtot;
  float Etotal;
  for (int i = 0; i < phoEta->size(); i++) {
    phoR = ecal_z/sinh(phoEta->at(i));
    phophi = phoPhi->at(i);
    esrhEtot = 0.0;
    hbherhEtot = 0.0;
    Etotal = 0.0;
    for (int j = 0; j < esrhE->size(); j++) {
      if (fabs(esrhEta->at(j)) < 1.566) continue;//only look at HE hits, not HB hits
      esrhR = sqrt(pow(esrhX->at(j),2)+pow(esrhY->at(j),2));
      esrhphi = esrhPhi->at(j);
      delphi = fabs(esrhphi-phophi);
      delr = fabs(phoR-esrhR);
      if (delphi > pi) delphi = 2*pi - delphi;
      if (delphi < 0.09 && (delr <= esrhZ->at(j)*tan(alpha))) esrhEtot+=esrhE->at(j); //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
    }
    for (int j = 0; j < hbherhE->size(); j++) {
      if (fabs(hbherhEta->at(j)) < 1.566) continue;//only look at HE hits, not HB hits
      hbherhR = sqrt(pow(hbherhX->at(j),2)+pow(hbherhY->at(j),2));
      hbherhphi = hbherhPhi->at(j);
      delphi = fabs(hbherhphi-phophi);
      delr = fabs(phoR-hbherhR);
      if (delphi > pi) delphi = 2*pi - delphi;
      if (delphi < 0.09 && (delr <= hbherhZ->at(j)*tan(alpha))) hbherhEtot+=hbherhE->at(j); //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
    }
    Etotal = esrhEtot + hbherhEtot;
    TotalVec.push_back(Etotal);
  }
  return TotalVec;
}

vector<float> phoHaloHE(vector<float> *esrhX, vector<float> *esrhY, vector<float> *esrhZ, vector<float> *esrhPhi, vector<float> *esrhE, vector<float> *esrhEta, vector<float> *hbherhX, vector<float> *hbherhY, vector<float> *hbherhZ, vector<float> *hbherhPhi, vector<float> *hbherhE, vector<float> *hbherhEta, vector<float> *phoEta, vector<float> *phoPhi) {
  vector<float> TotalVec;
  float pi = 3.14159265358979323846;
  float ecal_z = 319.5;//cm
  float alpha = pi/180;//1 degree in radians
  float phoR;
  float phophi;
  float phoeta;
  float phox;
  float phoy;
  float hbherhR;
  float hbherhphi;
  float hbherheta;
  float hbherhx;
  float hbherhy;
  //float delphi;
  //float delr;
  float dely;
  float delx;
  float del;
  float herhEtot;
  float Etotal;
  for (int i = 0; i < phoEta->size(); i++) {
    herhEtot = 0.0;
    phoeta = phoEta->at(i);
    if (fabs(phoeta)>=1.65 and fabs(phoeta)<=1.8) {//only encap photons get non-zero HaloHE
      phophi = phoPhi->at(i);
      phoR = ecal_z/sinh(phoEta->at(i));
      phox = phoR*cos(phophi);
      phoy = phoR*sin(phophi);
      for (int j = 0; j < hbherhE->size(); j++) {
        hbherheta = hbherhEta->at(j);
        if (fabs(hbherheta) < 1.65) continue;//only look at HE hits, not HB hits
        //hbherhR = sqrt(pow(hbherhX->at(j),2)+pow(hbherhY->at(j),2));
        //hbherhphi = hbherhPhi->at(j);
        //delphi = fabs(hbherhphi-phophi);
        hbherhx = hbherhX->at(j);
        hbherhy = hbherhY->at(j);
        delx = fabs(phox-hbherhx);
        dely = fabs(phoy-hbherhy);
        del = sqrt(pow(delx,2)+pow(dely,2));
        //if (delphi > pi) delphi = 2*pi - delphi;
        if (del <= 6 and check_samesign_eta(phoeta,hbherheta)) herhEtot+=hbherhE->at(j); //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
      }
    }
    TotalVec.push_back(herhEtot);
  }
  return TotalVec;
}

vector<float> phoHaloPre(vector<float> *esrhX, vector<float> *esrhY, vector<float> *esrhZ, vector<float> *esrhPhi, vector<float> *esrhE, vector<float> *esrhEta, vector<float> *hbherhX, vector<float> *hbherhY, vector<float> *hbherhZ, vector<float> *hbherhPhi, vector<float> *hbherhE, vector<float> *hbherhEta, vector<float> *phoEta, vector<float> *phoPhi) {
  vector<float> TotalVec;
  float pi = 3.14159265358979323846;
  float ecal_z = 319.5;//cm
  float alpha = pi/180;//1 degree in radians
  float phoR;
  float phophi;
  float phoeta;
  float phox;
  float phoy;
  float esrhR;
  float esrhphi;
  float esrheta;
  float esrhx;
  float esrhy;
  //float delphi;
  //float delr;
  float dely;
  float delx;
  float del;
  float esrhEtot;
  float Etotal;
  for (int i = 0; i < phoEta->size(); i++) {
    esrhEtot = 0.0;
    phoeta = phoEta->at(i);
    if (fabs(phoeta)>=1.566 and fabs(phoeta)<=2.6) {//only encap photons get non-zero HaloPre
      phophi = phoPhi->at(i);
      phoR = ecal_z/sinh(phoEta->at(i));
      phox = phoR*cos(phophi);
      phoy = phoR*sin(phophi);
      for (int j = 0; j < esrhE->size(); j++) {
        esrheta = esrhEta->at(j);
        if (fabs(esrheta) < 1.566) continue;//only look at HE hits, not HB hits
        //esrhR = sqrt(pow(esrhX->at(j),2)+pow(esrhY->at(j),2));
        //esrhphi = esrhPhi->at(j);
        //delphi = fabs(esrhphi-phophi);
        esrhx = esrhX->at(j);
        esrhy = esrhY->at(j);
        delx = fabs(phox-esrhx);
        dely = fabs(phoy-esrhy);
        del = sqrt(pow(delx,2)+pow(dely,2));
        //if (delphi > pi) delphi = 2*pi - delphi;
        if (del <= 3.1 and check_samesign_eta(phoeta,esrheta)) esrhEtot+=esrhE->at(j); //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
      }
    }
    TotalVec.push_back(esrhEtot);
  }
  return TotalVec;
}

bool near_object(float esrhx, float esrhy, float esrheta, vector<float>* obPhi, vector<float>* obEta, float ecal_z) {
  bool near = false;
  float obR;
  float obphi;
  float obeta;
  float obx;
  float oby;
  float dely;
  float delx;
  float del;
  int nob = obEta->size();
    for (int i = 0; i < nob; i++) {
      obphi = obPhi->at(i);
      obeta = obEta->at(i);
      obR = ecal_z/sinh(obeta);
      obx = obR*cos(obphi);
      oby = obR*sin(obphi);
      delx = fabs(obx-esrhx);
      dely = fabs(oby-esrhy);
      del = sqrt(pow(delx,2)+pow(dely,2));
      if (del <=13.5 and check_samesign_eta(obeta,esrheta) and fabs(obeta) > 1.65) near = true; // 13.5 cm is estimated from x,y plots of preshower hits and photons.
    }
  return near;
}

vector<int> esrhHaloPreID(vector<float> *esrhX, vector<float> *esrhY, vector<float> *esrhZ, vector<float> *esrhPhi, vector<float> *esrhE, vector<float> *esrhEta, vector<float> *hbherhX, vector<float> *hbherhY, vector<float> *hbherhZ, vector<float> *hbherhPhi, vector<float> *hbherhE, vector<float> *hbherhEta, vector<float> *phoEta, vector<float> *phoPhi, vector<float> *eleEta, vector<float> *elePhi, vector<float> *AK4CHSJet_Eta, vector<float> *AK4CHSJet_Phi, vector<float> *muEta, vector<float> *muPhi) {
  vector<int> esrhIDvec;
  float ecal_z = 319.5;//cm
  float esrhR;
  float esrhphi;
  float esrheta;
  float esrhx;
  float esrhy;
  int nesrh = esrhE->size();
  int id;
  for (int j = 0; j < nesrh; j++) {//esrh loop
    esrheta = esrhEta->at(j);
    esrhx = esrhX->at(j);
    esrhy = esrhY->at(j);
    id = 0;
    if (near_object(esrhx,esrhy,esrheta,phoPhi,phoEta,ecal_z)) id+=1;//near a photon
    if (near_object(esrhx,esrhy,esrheta,elePhi,eleEta,ecal_z)) id+=2;//near a electron
    if (near_object(esrhx,esrhy,esrheta,muPhi,muEta,ecal_z)) id+=4;//near a muon
    esrhIDvec.push_back(id);
  }
  return esrhIDvec;
}

int main(int argc, char** argv) {
  string infile = argv[1];
  TFile f(infile.c_str()); // Create TFile object for read file
  TTree* tree = (TTree*)f.Get("EventTree");
  tree->SetBranchStatus("*",1);

  string outname = barename(infile)+"_features.root";
  TFile fnew(outname.c_str(),"RECREATE"); // Write over last output file of same name if it exists
  auto newtree = tree->CloneTree(0);

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
  vector<float>* phoSCPhiWidth = nullptr;
  vector<float>* phoSCEtaWidth = nullptr;
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
  vector<float>* eleEta = nullptr;
  vector<float>* elePhi = nullptr;
  vector<float>* AK4CHSJet_Eta = nullptr;
  vector<float>* AK4CHSJet_Phi = nullptr;
  vector<float>* muEta = nullptr;
  vector<float>* muPhi = nullptr;
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
  tree->SetBranchAddress("eleEta",&eleEta);
  tree->SetBranchAddress("elePhi",&elePhi);
  tree->SetBranchAddress("AK4CHSJet_Eta",&AK4CHSJet_Eta);
  tree->SetBranchAddress("AK4CHSJet_Phi",&AK4CHSJet_Phi);
  tree->SetBranchAddress("muEta",&muEta);
  tree->SetBranchAddress("muPhi",&muPhi);
  vector<float> phoNumHERHnew;
  vector<float> phoNumHERHzsidenew;
  vector<float> phoNumESRHnew;
  vector<float> phoNumESRHzsidenew;
  vector<float> phoAHETotalnew;
  vector<float> phoHaloHEnew;
  vector<float> phoHaloPrenew;
  vector<int> esrhHaloPreIDnew;
  newtree->Branch("phoNumHERH",&phoNumHERHnew);
  newtree->Branch("phoNumHERHzside",&phoNumHERHzsidenew);
  newtree->Branch("phoNumESRH",&phoNumESRHnew);
  newtree->Branch("phoNumESRHzside",&phoNumESRHzsidenew);
  newtree->Branch("phoAHETotal",&phoAHETotalnew);
  newtree->Branch("phoHaloHE",&phoHaloHEnew);
  newtree->Branch("phoHaloPre",&phoHaloPrenew);
  newtree->Branch("esrhHaloPreID",&esrhHaloPreIDnew);
  int n = tree->GetEntries();
  float fn = n;
  for (Int_t i=0; i<n; i++) {
    tree->GetEntry(i);
    for (int j=0; j<phoEta->size(); j++) {
      phoNumHERHnew.push_back(get_close_phi_herh(phoPhi->at(j),hbherhPhi,hbherhEta));
      phoNumHERHzsidenew.push_back(get_close_phi_herhzside(phoPhi->at(j),phoEta->at(j),hbherhPhi,hbherhEta));
      phoNumESRHnew.push_back(get_close_phi_esrh(phoPhi->at(j),esrhPhi,esrhEta));
      phoNumESRHzsidenew.push_back(get_close_phi_esrhzside(phoPhi->at(j),phoEta->at(j),esrhPhi,esrhEta));
      phoAHETotalnew = phoAHETotal(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi);
      phoHaloHEnew = phoHaloHE(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi);
      phoHaloPrenew = phoHaloPre(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi);
      esrhHaloPreIDnew = esrhHaloPreID(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi, eleEta, elePhi, AK4CHSJet_Eta, AK4CHSJet_Phi, muEta, muPhi);
    }
    newtree->Fill();
    phoNumHERHnew.clear();
    phoNumHERHzsidenew.clear();
    phoNumESRHnew.clear();
    phoNumESRHzsidenew.clear();
    phoAHETotalnew.clear();
    phoHaloHEnew.clear();
    phoHaloPrenew.clear();
    esrhHaloPreIDnew.clear();
  }
  newtree->Write();
  fnew.Close();
  f.Close();
  return 0;
}
