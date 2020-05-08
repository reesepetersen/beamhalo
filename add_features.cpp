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

bool samesign(float thing1, float thing2) {
  if (thing1 > 0 && thing2 > 0) return true;
  else if (thing1 < 0 && thing2 < 0) return true;
  else return false;
}

//int get_close_phi_herh(float phoPhi,vector<float> *hbherhPhi,vector<float> *hbherhEta) {
//  int herh_count = 0;
//  float delphi;
//  float pi = 3.1415926535897932;
//  for (int i = 0; i < hbherhPhi->size(); i++) {
//    if (fabs(hbherhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
//    delphi = fabs(hbherhPhi->at(i)-phoPhi);
//    if (delphi > pi) delphi = 2*pi - delphi;
//    if (delphi < 0.09) herh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
//  }
//  return herh_count;
//}
//            
//int get_close_phi_herhzside(float phoPhi, float phoEta, vector<float> *hbherhPhi, vector<float> *hbherhEta) {
//  int herh_count = 0;
//  float delphi;
//  float pi = 3.1415926535897932;
//  for (int i = 0; i < hbherhPhi->size(); i++) {
//    if (fabs(hbherhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
//    delphi = fabs(hbherhPhi->at(i)-phoPhi);
//    if (delphi > pi) delphi = 2*pi - delphi;
//    bool samesign_eta = samesign(phoEta,hbherhEta->at(i));
//    if (delphi < 0.09 && samesign_eta) herh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
//  }
//  return herh_count;
//}
//            
//int get_close_phi_esrh(float phoPhi,vector<float> *esrhPhi,vector<float> *esrhEta) {
//  int esrh_count = 0;
//  float delphi;
//  float pi = 3.1415926535897932;
//  for (int i = 0; i < esrhPhi->size(); i++) {
//    if (fabs(esrhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
//    delphi = fabs(esrhPhi->at(i)-phoPhi);
//    if (delphi > pi) delphi = 2*pi - delphi;
//    if (delphi < 0.09) esrh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
//  }
//  return esrh_count;
//}
//
//int get_close_phi_esrhzside(float phoPhi, float phoEta, vector<float> *esrhPhi,vector<float> *esrhEta) {
//  int esrh_count = 0;
//  float delphi;
//  float pi = 3.1415926535897932;
//  for (int i = 0; i < esrhPhi->size(); i++) {
//    if (fabs(esrhEta->at(i)) < 1.566) continue;//only look at HE hits, not HB hits
//    delphi = fabs(esrhPhi->at(i)-phoPhi);
//    if (delphi > pi) delphi = 2*pi - delphi;
//    bool samesign_eta = samesign(phoEta,esrhEta->at(i));
//    if (delphi < 0.09 && samesign_eta) esrh_count++; //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
//  }
//  return esrh_count;
//}

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
        if (del <= 25 and samesign(phoeta,hbherheta)) herhEtot+=hbherhE->at(j); //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
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
        if (del <= 13.5 and samesign(phoeta,esrheta)) esrhEtot+=esrhE->at(j); //0.09 radians is just over 5 degrees. Most HE segements cover 10 degrees in phi
      }
    }
    TotalVec.push_back(esrhEtot);
  }
  return TotalVec;
}

bool near_object(float rhx, float rhy, float rheta, vector<float>* obPhi, vector<float>* obEta, float ecal_z, float rad) {
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
      delx = fabs(obx-rhx);
      dely = fabs(oby-rhy);
      del = sqrt(pow(delx,2)+pow(dely,2));
      if (del <=rad and samesign(obeta,rheta) and fabs(obeta) > 1.65) near = true;
    }
  return near;
}

//vector<int> esrhHaloPreID(vector<float> *esrhX, vector<float> *esrhY, vector<float> *esrhZ, vector<float> *esrhPhi, vector<float> *esrhE, vector<float> *esrhEta, vector<float> *hbherhX, vector<float> *hbherhY, vector<float> *hbherhZ, vector<float> *hbherhPhi, vector<float> *hbherhE, vector<float> *hbherhEta, vector<float> *phoEta, vector<float> *phoPhi, vector<float> *eleEta, vector<float> *elePhi, vector<float> *AK4CHSJet_Eta, vector<float> *AK4CHSJet_Phi, vector<float> *muEta, vector<float> *muPhi) {
//  vector<int> esrhIDvec;
//  float ecal_z = 319.5;//cm
//  float rad = 13.5;//cm
//  float esrhR;
//  float esrhphi;
//  float esrheta;
//  float esrhx;
//  float esrhy;
//  int nesrh = esrhE->size();
//  int id;
//  for (int j = 0; j < nesrh; j++) {//esrh loop
//    esrheta = esrhEta->at(j);
//    esrhx = esrhX->at(j);
//    esrhy = esrhY->at(j);
//    id = 0;
//    if (near_object(esrhx,esrhy,esrheta,phoPhi,phoEta,ecal_z,rad)) id+=1;//near a photon
//    if (near_object(esrhx,esrhy,esrheta,elePhi,eleEta,ecal_z,rad)) id+=2;//near a electron
//    if (near_object(esrhx,esrhy,esrheta,muPhi,muEta,ecal_z,rad)) id+=4;//near a muon
//    esrhIDvec.push_back(id);
//  }
//  return esrhIDvec;
//}
//
//vector<int> herhHaloHEID(vector<float> *esrhX, vector<float> *esrhY, vector<float> *esrhZ, vector<float> *esrhPhi, vector<float> *esrhE, vector<float> *esrhEta, vector<float> *hbherhX, vector<float> *hbherhY, vector<float> *hbherhZ, vector<float> *hbherhPhi, vector<float> *hbherhE, vector<float> *hbherhEta, vector<float> *phoEta, vector<float> *phoPhi, vector<float> *eleEta, vector<float> *elePhi, vector<float> *AK4CHSJet_Eta, vector<float> *AK4CHSJet_Phi, vector<float> *muEta, vector<float> *muPhi) {
//  vector<int> hbherhIDvec;
//  float ecal_z = 319.5;//cm
//  float rad = 25;//cm
//  float hbherhR;
//  float hbherhphi;
//  float hbherheta;
//  float hbherhx;
//  float hbherhy;
//  int nhbherh = hbherhE->size();
//  int id;
//  for (int j = 0; j < nhbherh; j++) {//herh loop
//    hbherheta = hbherhEta->at(j);
//    hbherhx = hbherhX->at(j);
//    hbherhy = hbherhY->at(j);
//    id = 0;
//    if (near_object(hbherhx,hbherhy,hbherheta,phoPhi,phoEta,ecal_z,rad)) id+=1;//near a photon
//    if (near_object(hbherhx,hbherhy,hbherheta,elePhi,eleEta,ecal_z,rad)) id+=2;//near a electron
//    if (near_object(hbherhx,hbherhy,hbherheta,muPhi,muEta,ecal_z,rad)) id+=4;//near a muon
//    hbherhIDvec.push_back(id);
//  }
//  return hbherhIDvec;
//}

vector<int> hbherhReg(float phoSCPhi, float phoSCEta, vector<float>* hbherhE, vector<float>* hbherhX, vector<float>* hbherhY, vector<float>* hbherhZ, vector<float>* hbherhPhi, vector<float>* hbherhEta) {
  float pi = 3.14159265358979323846;
  float ecal_z = 319.5;
  int nhbherh = hbherhE->size();
  float dphi;
  float dx;
  float dy;
  float phoSCR = ecal_z/fabs(sinh(phoSCEta));
  float phox = phoSCR*cos(phoSCPhi);
  float phoy = phoSCR*sin(phoSCPhi);
  float drho;
  float deta;
  float dR;
  int hbherhR;
  vector<int> Reg;
  for(int i = 0; i < nhbherh; i++) {
    hbherhR = sqrt(pow(hbherhX->at(i),2)+pow(hbherhY->at(i),2));
    if (not samesign(phoSCEta,hbherhEta->at(i)) or fabs(phoSCEta)<1.65 or fabs(phoSCEta)>1.8) {//eta restriction to be removed after testing
      Reg.push_back(0);
      continue;
    }
    dphi = fabs(hbherhPhi->at(i)-phoSCPhi);
    if (dphi > pi) dphi = 2*pi - dphi;
    dx = fabs(hbherhX->at(i)-phox);
    dy = fabs(hbherhY->at(i)-phoy);
    drho = sqrt(pow(dx,2)+pow(dy,2));
    deta = fabs(hbherhEta->at(i)-phoSCEta);
    dR = sqrt(pow(dphi,2)+pow(deta,2));
    if (dphi < 0.15 and drho < 26) Reg.push_back(2);
    else if (dR <= 0.15) Reg.push_back(1);
    else Reg.push_back(0);
  }
  return Reg;
}

vector<float> E1A1E2A2(float phoSCPhi, float phoSCEta, vector<float>* hbherhE, vector<float>* hbherhX, vector<float>* hbherhY, vector<float>* hbherhZ, vector<int> Reg) {
  int reg_check = 0;
  for (int i = 0; i < Reg.size(); i++) {
    reg_check += Reg[i];
  }
  if (reg_check == 0 or fabs(phoSCEta) < 1.65 or fabs(phoSCEta) > 1.8) return {-4,-4,-4,-4}; //temporary, but keep this restriction for barrel photons.
  vector<float> EA1EA2;
  float ecal_z = 319.5;
  float phoSCR = ecal_z/fabs(sinh(phoSCEta));
  float phoX = phoSCR*cos(phoSCPhi);
  float phoY = phoSCR*sin(phoSCPhi);
  float phoZ = ecal_z;
  int nhbherh = hbherhE->size();
  float E;
  float E1 = 0;
  float E2 = 0;
  float EX1 = 0;
  float EX2 = 0;
  float EY1 = 0;
  float EY2 = 0;
  float EZ1 = 0;
  float EZ2 = 0;
  float A1;
  float A2;
  int reg;
  for(int i = 0; i < nhbherh; i++) {
    reg = Reg[i];
    if (reg == 0) continue;
    if (reg == 1) {
      E = hbherhE->at(i);
      E1 += E;
      EX1 += E*hbherhX->at(i);
      EY1 += E*hbherhY->at(i);
      EZ1 += E*hbherhZ->at(i);
    }
    else if (reg == 2) {
      E = hbherhE->at(i);
      E2 += E;
      EX2 += E*hbherhX->at(i);
      EY2 += E*hbherhY->at(i);
      EZ2 += E*hbherhZ->at(i);
    }
  }
  if (E1 > 0) {
    float avx1 = EX1/E1;
    float avy1 = EY1/E1;
    float avr1 = sqrt(pow(avx1,2)+pow(avy1,2));
    float dr1 = avr1-phoSCR;
    float avz1 = EZ1/E1;
    float dz1 = fabs(avz1)-phoZ;
    A1 = atan2(dr1,dz1);
  }
  else {
    A1 = -4;
  }
  if (E2 > 0) {
    float avx2 = EX2/E2;
    float avy2 = EY2/E2;
    float avr2 = sqrt(pow(avx2,2)+pow(avy2,2));
    float avz2 = EZ2/E2;
    float dr2 = avr2-phoSCR;
    float dz2 = fabs(avz2-phoZ);
    A2 = atan2(dr2,dz2);
  }
  else {
    A2 = -4;
  }
  EA1EA2.push_back(E1);
  EA1EA2.push_back(A1);
  EA1EA2.push_back(E2);
  EA1EA2.push_back(A2);
  return EA1EA2;
}

float phoEnWeAn(vector<float> EA12) {
  float E1 = EA12[0];
  float A1 = EA12[1];
  float E2 = EA12[2];
  float A2 = EA12[3];
  if (E1 < 0) return 4;
  float enwean = (A1*E1+A2*E2)/(E1+E2);
  return enwean;
}

vector<int> esrhReg(float phoSCPhi, float phoSCEta, vector<float>* esrhE, vector<float>* esrhX, vector<float>* esrhY, vector<float>* esrhZ, vector<float>* esrhPhi, vector<float>* esrhEta, vector<float>* esrhPlane) {
  float pi = 3.14159265358979323846;
  float ecal_z = 319.5;
  int nesrh = esrhE->size();
  float dx;
  float dy;
  float dfx1;
  float dfy1;
  float dfx2;
  float dfy2;
  float phoSCR = ecal_z/fabs(sinh(phoSCEta));
  float phox = phoSCR*cos(phoSCPhi);
  float phoy = phoSCR*sin(phoSCPhi);
  float phoTheta = 2*atan(exp(-phoSCEta));
  float phoZAngle;
  if (phoTheta > pi/2) {
    phoZAngle = pi - phoTheta;
  }
  else {
    phoZAngle = phoTheta;
  }
  float halfAngle = phoZAngle/2;
  float fauSCR1 = phoSCR-15.3/sinh(fabs(phoSCEta));
  float fauSCR2 = phoSCR-10.65/sinh(fabs(phoSCEta));
  float halfanR1 = phoSCR-15.3*tan(halfAngle)/2;
  float halfanR2 = phoSCR-10.65*tan(halfAngle)/2;
  //float faux1 = fauSCR1*cos(phoSCPhi);
  //float fauy1 = fauSCR1*sin(phoSCPhi);
  //float faux2 = fauSCR2*cos(phoSCPhi);
  //float fauy2 = fauSCR2*sin(phoSCPhi);
  float drho;
  float drho1;
  float drho2;
  float esrhR;
  float dphi;
  float deta;
  float dR;
  int reg;
  vector<int> Reg;
  for(int i = 0; i < nesrh; i++) {
    if (fabs(phoSCEta) < 1.65 or fabs(phoSCEta) > 1.8 or not samesign(phoSCEta,esrhZ->at(i))) {
      Reg.push_back(0);
      continue;
    }
    esrhR = sqrt(pow(esrhX->at(i),2)+pow(esrhY->at(i),2));
    dx = fabs(esrhX->at(i)-phox);
    dy = fabs(esrhY->at(i)-phoy);
    //dfx1 = fabs(esrhX->at(i)-faux1);
    //dfy1 = fabs(esrhY->at(i)-fauy1);
    //dfx2 = fabs(esrhX->at(i)-faux2);
    //dfy2 = fabs(esrhY->at(i)-fauy2);
    drho = sqrt(pow(dx,2)+pow(dy,2));
    //drho1 = sqrt(pow(dfx1,2)+pow(dfy1,2));
    //drho2 = sqrt(pow(dfx2,2)+pow(dfy2,2));
    dphi = fabs(esrhPhi->at(i)-phoSCPhi);
    deta = fabs(esrhEta->at(i)-phoSCEta);
    dR = sqrt(pow(dphi,2)+pow(deta,2));
    reg = 0;
    if (drho < 2 and ((esrhR >= halfanR1 and esrhPlane->at(i) == 1) or (esrhR >= halfanR2 and esrhPlane->at(i) == 2))) reg = 2;
    //if ((drho1 <= 13.5 and esrhR < halfanR1 and esrhPlane->at(i)==1) or (drho2 <= 13.5 and esrhR < halfanR2 and esrhPlane->at(i)==2)) reg += 1;
    if (dR < 0.017 and ((esrhR < halfanR1 and esrhPlane->at(i)==1) or (esrhR < halfanR2 and esrhPlane->at(i)==2))) reg = 1;
    Reg.push_back(reg);
  }
  return Reg;
}

vector<float> preE1A1E2A2(float phoSCPhi, float phoSCEta, vector<float>* esrhE, vector<float>* esrhX, vector<float>* esrhY, vector<float>* esrhZ, vector<int> Reg) {
  float pi = 3.14159265358979323846;
  int reg_check = 0;
  for (int i = 0; i < Reg.size(); i++) {
    reg_check += Reg[i];
  }
  if (reg_check == 0 or fabs(phoSCEta) < 1.65 or fabs(phoSCEta) > 1.8) return {-4,-4,-4,-4}; //temporary, but keep this restriction for barrel photons.
  vector<float> EA1EA2;
  float ecal_z = 319.5;
  float phoSCR = ecal_z/fabs(sinh(phoSCEta));
  float phoX = phoSCR*cos(phoSCPhi);
  float phoY = phoSCR*sin(phoSCPhi);
  float phoZ = ecal_z;
  float phoTheta = 2*atan(exp(-phoSCEta));
  float phoZAngle;
  if (phoTheta > pi/2) {
    phoZAngle = pi - phoTheta;
  }
  float halfAngle = phoZAngle/2;
  int nesrh = esrhE->size();
  float E;
  float E1 = 0;
  float E2 = 0;
  float EX1 = 0;
  float EX2 = 0;
  float EY1 = 0;
  float EY2 = 0;
  float EZ1 = 0;
  float EZ2 = 0;
  float A1;
  float A2;
  int reg;
  for(int i = 0; i < nesrh; i++) {
    reg = Reg[i];
    if (reg == 0) continue;
    if (reg == 1) {
      E = esrhE->at(i);
      E1 += E;
      EX1 += E*esrhX->at(i);
      EY1 += E*esrhY->at(i);
      EZ1 += E*esrhZ->at(i);
    }
    if (reg == 2) {
      E = esrhE->at(i);
      E2 += E;
      EX2 += E*esrhX->at(i);
      EY2 += E*esrhY->at(i);
      EZ2 += E*esrhZ->at(i);
    }
  }
  if (E1 > 0) {
    float avx1 = EX1/E1;
    float avy1 = EY1/E1;
    float avr1 = sqrt(pow(avx1,2)+pow(avy1,2));
    float dr1 = phoSCR-avr1;
    float avz1 = EZ1/E1;
    float dz1 = phoZ-fabs(avz1);
    A1 = atan2(dr1,dz1);
  }
  else {
    A1 = -4;
  }
  if (E2 > 0) {
    float avx2 = EX2/E2;
    float avy2 = EY2/E2;
    float avr2 = sqrt(pow(avx2,2)+pow(avy2,2));
    float dr2 = phoSCR-avr2;
    float avz2 = EZ2/E2;
    float dz2 = phoZ-fabs(avz2);
    A2 = atan2(dr2,dz2);
  }
  else {
    A2 = -4;
  }
  EA1EA2.push_back(E1);
  EA1EA2.push_back(A1);
  EA1EA2.push_back(E2);
  EA1EA2.push_back(A2);
  return EA1EA2;
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
  //vector<float> phoNumHERHnew;
  //vector<float> phoNumHERHzsidenew;
  //vector<float> phoNumESRHnew;
  //vector<float> phoNumESRHzsidenew;
  vector<float> phoAHETotalnew;
  vector<float> phoHaloHEnew;
  vector<float> phoHaloPrenew;
  vector<float> phoE1new;
  vector<float> phoA1new;
  vector<float> phoE2new;
  vector<float> phoA2new;
  vector<float> phoEnWeAnnew;
  vector<float> phoPreE1new;
  vector<float> phoPreA1new;
  vector<float> phoPreE2new;
  vector<float> phoPreA2new;
  vector<float> phoPreEnWeAnnew;
  //vector<int> esrhHaloPreIDnew;
  //vector<int> herhHaloHEIDnew;
  newtree->Branch("phoAHETotal",&phoAHETotalnew);
  newtree->Branch("phoHaloHE",&phoHaloHEnew);
  newtree->Branch("phoHaloPre",&phoHaloPrenew);
  newtree->Branch("phoE1",&phoE1new);
  newtree->Branch("phoA1",&phoA1new);
  newtree->Branch("phoE2",&phoE2new);
  newtree->Branch("phoA2",&phoA2new);
  newtree->Branch("phoEnWeAn",&phoEnWeAnnew);
  newtree->Branch("phoPreE1",&phoPreE1new);
  newtree->Branch("phoPreA1",&phoPreA1new);
  newtree->Branch("phoPreE2",&phoPreE2new);
  newtree->Branch("phoPreA2",&phoPreA2new);
  newtree->Branch("phoPreEnWeAn",&phoPreEnWeAnnew);
  //newtree->Branch("esrhHaloPreID",&esrhHaloPreIDnew);
  //newtree->Branch("herhHaloHEID",&herhHaloHEIDnew);
  int n = tree->GetEntries();
  float fn = n;
  for (Int_t i=0; i<n; i++) {
    tree->GetEntry(i);
    for (int j=0; j<phoEta->size(); j++) {
      vector<int> heReg = hbherhReg(phoSCPhi->at(j), phoSCEta->at(j), hbherhE, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhEta);
      vector<int> preReg = esrhReg(phoSCPhi->at(j), phoSCEta->at(j), esrhE, esrhX, esrhY, esrhZ, esrhPhi, esrhEta, esrhPlane);
      vector<float> EA12 = E1A1E2A2(phoSCPhi->at(j), phoSCEta->at(j), hbherhE, hbherhX, hbherhY, hbherhZ, heReg);
      vector<float> preEA12 = preE1A1E2A2(phoSCPhi->at(j), phoSCEta->at(j), esrhE, esrhX, esrhY, esrhZ, preReg);
      phoE1new.push_back(EA12[0]);
      phoA1new.push_back(EA12[1]);
      phoE2new.push_back(EA12[2]);
      phoA2new.push_back(EA12[3]);
      phoEnWeAnnew.push_back(phoEnWeAn(EA12));
      phoPreE1new.push_back(preEA12[0]);
      phoPreA1new.push_back(preEA12[1]);
      phoPreE2new.push_back(preEA12[2]);
      phoPreA2new.push_back(preEA12[3]);
      phoPreEnWeAnnew.push_back(phoEnWeAn(preEA12));
    }
    phoAHETotalnew = phoAHETotal(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi);
    phoHaloHEnew = phoHaloHE(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi);
    phoHaloPrenew = phoHaloPre(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi);
    //esrhHaloPreIDnew = esrhHaloPreID(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi, eleEta, elePhi, AK4CHSJet_Eta, AK4CHSJet_Phi, muEta, muPhi);
    //herhHaloHEIDnew = herhHaloHEID(esrhX, esrhY, esrhZ, esrhPhi, esrhE, esrhEta, hbherhX, hbherhY, hbherhZ, hbherhPhi, hbherhE, hbherhEta, phoEta, phoPhi, eleEta, elePhi, AK4CHSJet_Eta, AK4CHSJet_Phi, muEta, muPhi);
    newtree->Fill();
    phoAHETotalnew.clear();
    phoHaloHEnew.clear();
    phoHaloPrenew.clear();
    phoE1new.clear();
    phoA1new.clear();
    phoE2new.clear();
    phoA2new.clear();
    phoEnWeAnnew.clear();
    phoPreE1new.clear();
    phoPreA1new.clear();
    phoPreE2new.clear();
    phoPreA2new.clear();
    phoPreEnWeAnnew.clear();
    //esrhHaloPreIDnew.clear();
    //herhHaloHEIDnew.clear();
  }
  newtree->Write();
  fnew.Close();
  f.Close();
  return 0;
}
