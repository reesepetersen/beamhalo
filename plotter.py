#!/usr/bin/env python

# Author : Reese Petersen
# Date   : 10 Dec 2019
# Inst   : University of Minnesota, CMS Collaboration

# import useful libraries
import os
import math
import sys
import argparse
import subprocess as sp
import ROOT as r
import datetime
import re
from collections import OrderedDict

# taking an argument from the command line
#parser = argparse.ArgumentParser(description = 'This runs a monte carlo simulation for the energy flux of cosmic muons at sea level on Earth as a function of energy.')
#parser.add_argument("--inFile" , dest="inFile" , help="File with ttree", required=True)
#arg = parser.parse_args()

print "Welcome to RapidPlot (TM)!\n"

rtypedict = {
    'C':'cstring',
    'B':'Char_t',
    'b':'UChar_t',
    'S':'Short_t',
    's':'UShort_t',
    'I':'Int_t',
    'i':'UInt_t',
    'F':'Float_t',
    'f':'Float16_t',
    'D':'Double_t',
    'd':'Double32_t',
    'L':'Long64_t',
    'l':'ULong64_t',
    'O':'Bool_t'}

# check input file for .root extension
#filestring = arg.inFile
filestring = sys.argv[1]
if filestring.split('.')[-1] != 'root':
  print "Input file does not have .root extentsion. exiting..."
  quit()

# open file and check if opened
#rf = r.TFile(filestring)
#if not rf.IsOpen():
  #print "Could not open file, exiting..."
  #quit()

def rootls(filepath,mode='-l'):
  out = sp.Popen(['rootls',mode,filepath],stdout=sp.PIPE,stderr=sp.STDOUT)
  stdout,stderr = out.communicate()
  if stderr != None:
    print stderr
    quit()
  return stdout

def noempty(instring):
  if instring != '':
    return True
  else:
    return False

def clean_list(optlist):
  cleanlist = []
  for opt in optlist:
    if opt not in cleanlist and opt != {}:
      cleanlist.append(opt)
  return cleanlist

def ls_to_list(rls):
  rlslist = rls.split('\n')
  return filter(noempty,rlslist)

def parse_ls(rls,mode='-l'):
  filelist = ls_to_list(rls)
  optlist = []
  for entry in filelist:
    entrylist = filter(noempty,re.split('[\s]{2,}',entry))
    entrydict = {'cpptype':entrylist[0],'name':entrylist[2]}
    optlist.append(entrydict)
  optlist = clean_list(optlist)
  return optlist

def parse(dirstring,mode='-l'):
  rlsdir = rootls(dirstring,mode)
  print dirstring+':'
  dirlist = parse_ls(rlsdir,mode)
  for item in dirlist:
    item['dir'] = dirstring
  return dirlist

def process_directory(dirstring):
  dirlist = parse(dirstring)
  chosen = {}
  if len(dirlist) > 1:
    for i in range(len(dirlist)):
      entry = dirlist[i]
      ctype = entry['cpptype']
      name = entry['name']
      print str(i)+': '+ctype+' '+name
    choice = raw_input("Please enter a number: ")
    chosen = dirlist[int(choice)]
  elif len(dirlist) == 1:
    chosen = dirlist[0]
    print chosen['cpptype']+' '+chosen['name']
  return chosen

def get_filestring(dirstring):
  if ':' in dirstring:
    return dirstring.split(':')[0]
  else:
    return dirstring

def get_barename(filestring):
  if '/' in filestring:
    return filestring.split('/')[-1].split('.')[0]
  else:
    return filestring.split('.')[0]

def get_treepath(dirstring,tree,joiner = '/'):
  if ':' in dirstring:
    dirlist = dirstring.split(':')
    dirlist.pop(0)
    return joiner.join(dirlist)+joiner+tree
  else:
    return tree

def parse_tree_entry(entrylist):
  if len(entrylist) == 3:
    nametype = entrylist[1]
    if '/' in nametype:
      t = nametype.split('"')[1].split('/')[1]
      rtype = rtypedict[t]
      return {'name':entrylist[0],'cpptype':rtype}
    elif '"'+entrylist[0]+'"' == entrylist[1]:
      return {'name':entrylist[0],'cpptype':'unknown'}
  else:
    return {}

def parse_ls_ttree(rls,dirstring,treename):# Only take in the entries after treename is found
  treelist = clean_list(ls_to_list(rls))
  optlist = []
  right_tree = False
  for entry in treelist:
    if "TTree" in entry:
      if re.split('[\s]{2,}',entry)[2] == treename:
        right_tree = True
      else:
        right_tree = False
    if right_tree:
      entrylist = filter(noempty,re.split('[\s]{2,}',entry))
      entrydict = parse_tree_entry(entrylist)
      optlist.append(entrydict)
  optlist = clean_list(optlist)
  unknown_detected = False
  for entry in optlist:
    if entry['cpptype'] == 'unknown':
      unknown_detected = True
      break
  if unknown_detected:
    rootfile = r.TFile(get_filestring(dirstring))
    tree = rootfile.Get(get_treepath(dirstring,treename))
    for entry in optlist:
      if entry['cpptype'] == 'unknown':
        entry['cpptype'] = tree.GetLeaf(entry['name']).GetTypeName()
    rootfile.Close()
  return optlist

def parse_ttree(dirstring,treename):# This function should only return the list of branches corresponding to treename
  rls = rootls(dirstring,'-t')
  treelist = parse_ls_ttree(rls,dirstring,treename)
  return treelist

def select_from_ttree(dirstring,treename):
  print "TTree: "+treename
  treelist = parse_ttree(dirstring,treename)
  chosen = {}
  for i in range(len(treelist)):
    entry = treelist[i]
    ctype = entry['cpptype']
    name = entry['name']
    print str(i)+': '+ctype+' '+name
  choice = raw_input("Please enter a number: ")
  chosen = treelist[int(choice)]
  print 'selected: '+chosen['cpptype']+' '+chosen['name']
  return chosen

def plot_treebranch(dirstring,treename):
  choice = select_from_ttree(dirstring,treename)
  rootfile = r.TFile(get_filestring(dirstring))
  tree = rootfile.Get(get_treepath(dirstring,treename))
  sel = raw_input("Selection (default none): ")
  if choice['name'] == 'phoNumHERH':
    hherh = r.TH1D('hherh','phoNumHERH',100,0,30)
    tree.Draw(choice['name']+'>>hherh',sel)
  else:
    tree.Draw(choice['name'],sel)
  raw_input("Enter to quit ")

def compare_var_from_files(name,dirstring,treename,compfile,hist_range,hist_units):
  hrange = hist_range.split(',')
  hrlow = float(hrange[0])
  hrhigh = float(hrange[1])
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf1 = r.TFile(get_filestring(dirstring))
  rf2 = r.TFile(compfile)
  t1 = rf1.Get(get_treepath(dirstring,treename))
  t2 = rf2.Get(get_treepath(dirstring,treename))
  h1 = r.TH1D('h1',name,100,hrlow,hrhigh)
  h2 = r.TH1D('h2',name,100,hrlow,hrhigh)
  if 'prompt' in dirstring:
    h1.SetLineColor(r.kRed)
    h2.SetLineColor(r.kBlue)
  else:
    h1.SetLineColor(r.kBlue)
    h2.SetLineColor(r.kRed)
  h1.GetXaxis().SetTitle(hist_units)
  sel = raw_input("Selection (default none): ")
  t1.Draw(name+'>>h1',sel,'goff')
  t2.Draw(name+'>>h2',sel,'goff')
  comp_canvas = r.TCanvas('comp_canvas','compare '+name)
  if name == "esrhE":
    comp_canvas.SetLogy()
  sum1 = h1.Integral()
  sum2 = h2.Integral()
  h1.Scale(1/sum1)
  h2.Scale(1/sum2)
  h1.Draw('hist')
  h2.Draw('same hist')
  sum1_text = r.TPaveText(0.1,0.9,0.5,1.0)
  sum1_text.AddText(get_filestring(dirstring))
  sum2_text = r.TPaveText(0.1,0.8,0.5,0.9)
  sum2_text.AddText(compfile)
  sum1_text.Draw();
  sum2_text.Draw();
  leg = r.TLegend(0.5,0.9,0.9,1.0)
  leg.AddEntry(h1,get_barename(get_filestring(dirstring)))
  leg.AddEntry(h2,get_barename(compfile))
  leg.Draw()
  comp_canvas.Update()
  saveyn = raw_input("Save cavnas? (y/n): ")
  if saveyn == 'y':
    comp_canvas.SaveAs(name+'_'+get_barename(get_filestring(dirstring))+'_'+get_barename(compfile)+'.png')

def plot_2D(x_var,x_range,x_units,y_var,y_range,y_units,dirstring,treename):
  xr = x_range.split(',')
  yr = y_range.split(',')
  polar = raw_input("Polar Plot? (y/n): ")
  if polar == 'y':
    xlow = -float(yr[1])
    ylow = xlow
    xhigh = float(yr[1])
    yhigh = xhigh
  else:
    xlow = float(xr[0])
    ylow = float(yr[0])
    xhigh = float(xr[1])
    yhigh = float(yr[1])
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_treepath(dirstring,treename))
  h = r.TH2F('h',x_var+" vs "+y_var,200,xlow,xhigh,200,ylow,yhigh)
  if polar != 'y':
    h.GetXaxis().SetTitle(x_units)
    h.GetYaxis().SetTitle(y_units)
  sel = raw_input("Selection (continue for none): ")
  if polar == 'y':
    t.Draw(y_var+":"+x_var+">>h",sel,'goff pol')
  else:
    t.Draw(y_var+":"+x_var+">>h",sel,'goff')
  corfac = h.GetCorrelationFactor()
  print "Pearson Correlation Factor: "+str(corfac)
  canv = r.TCanvas('canvas2D','2D '+x_var+' vs '+y_var)
  if polar == 'y':
    h.Draw('colz pol')
  else:
    h.Draw('colz')
  canv.Update()
  save = raw_input("Save? (y/n): ")
  if save == 'y':
    canv.SaveAs(x_var+'_vs_'+y_var+'_'+get_barename(get_filestring(dirstring))+'.png')

def plot_he_r(dirstring,treename):
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_treepath(dirstring,treename))
  t.Draw("sqrt(hbherhX*hbherhX+hbherhY*hbherhY)","hbherhEta>1.48 && hbherhE>1.5")
  raw_input('wait')

def plot_HE_phoSC_prox(dirstring,treename):
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_treepath(dirstring,treename))
  t.Draw("sqrt((hbherhX-phoSCX)^2+(hbherhY-phoSCY)^2)","hbherhE>1 && hbherhEta>1.3")
  raw_input('wait')

def compare_var(dirstring,treename):
  choice = select_from_ttree(dirstring,treename)
  compfile = raw_input("Please give second file for comparison: ")
  hist_range = raw_input("Please give the range (e.g. -3.14,3.14): ")
  hist_units = raw_input("Please give the units: ")
  compare_var_from_files(choice['name'],dirstring,treename,compfile,hist_range,hist_units)

def process_x_v_y(dirstring,treename):
  choice = select_from_ttree(dirstring,treename)
  x_var = choice['name']
  print "x variable: "+x_var
  x_range = raw_input("Please give the x range (e.g. -3.14,3.14): ")
  x_units = raw_input("Please give the x units: ")
  y_var = raw_input("y variable: ")
  y_range = raw_input("Please give the y range (e.g. -3.14,3.14): ")
  y_units = raw_input("Please give the y units: ")
  plot_2D(x_var,x_range,x_units,y_var,y_range,y_units,dirstring,treename)

def make_prompt_filter(dirstring,treename):
  selector_file = open("prompt_filter.cpp",'w')
  selector_file.write('''#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <numeric>
#include "TFile.h"
#include "TTree.h"
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

int max_index (vector<float> vin) {
  int maxi = -1;
  int max = 0;
  int num = vin.size();
  if (num > 0) {
    max = vin.at(0);
    maxi = 0;
    for (int i=1; i<num; i++) {
      if (vin.at(i) > max) {
        max = vin.at(i);
        maxi = i;
      }
    }
  }
  return maxi;
}

main(int argc, char** argv) {
  string inputfilename = argv[1];
  ifstream inputfile;
  string inputcontent;
  inputfile.open(inputfilename);
  getline(inputfile,inputcontent);
  inputfile.close();
  string infile = inputcontent;
  TFile f(infile.c_str());
  TTree* tree = (TTree*)f.Get("EventTree");
  tree->SetBranchStatus("*",1);
  string outname = barename(infile)+"_promptcut.root";
  TFile fnew(outname.c_str(),"RECREATE");
  TTree newtree("EventTree","Event data");
''')
  treelist = parse_ttree(dirstring,treename)
  ctypes = []
  names = []
  for i in range(len(treelist)):
    entry = treelist[i]
    ctypes.append(entry['cpptype'])
    names.append(entry['name'])
  for i in range(len(treelist)):# infile branch initialization
    if 'vector' in ctypes[i]:
      selector_file.write('  '+ctypes[i]+'* '+names[i]+' = nullptr;\n')
    else:
      selector_file.write('  '+ctypes[i]+' '+names[i]+';\n')
  for i in range(len(treelist)):# infile SetBranchAddress
    selector_file.write('  tree->SetBranchAddress("'+names[i]+'",&'+names[i]+');\n')
  for i in range(len(treelist)):# outfile branch initialization
    selector_file.write('  '+ctypes[i]+' '+names[i]+'new;\n')
  for i in range(len(treelist)):# outfile branch initialization
    selector_file.write('  newtree.Branch("'+names[i]+'",&'+names[i]+'new);\n')
  selector_file.write('''  float prompt_pass_count = 0.0f;
  float pho_total_count = 0.0f;
  int n = tree->GetEntries();
  float fn = n;
  float phoEta_pass_count = 0.0f;
  float phoHoverE_pass_count = 0.0f;
  float phoPhi_pass_count = 0.0f;
  float phoPFChWorstIso_pass_count = 0.0f;
  float phoEt_pass_count = 0.0f;
  bool first_time;
  for (Int_t i=0; i<n; i++) {
    tree->GetEntry(i);
    first_time = true;
    for (int j=0; j<phoEta->size(); j++) {
      if(fabs(phoEta->at(j)) > 1.566 && fabs(phoEta->at(j)) < 2.5) {
        if(phoHoverE->at(j) < 0.026) {
          if(phoPFChWorstIso->at(j) < 1.146) {
            if (fabs(phoPhi->at(j)) < 2.8 && fabs(phoPhi->at(j)) > 0.2 ) {
              if (phoEt->at(j) > 250.0) {
''')
  for i in range(len(treelist)):
    if 'pho' in names[i] and 'vector' in ctypes[i]:# keep photon info only for passing photons
      selector_file.write('                '+names[i]+'new.push_back('+names[i]+'->at(j));\n')
  selector_file.write('                if (first_time) {\n')
  for i in range(len(treelist)):# the first time a passing photon is found in an event, save all non-pho vectors
    if 'pho' not in names[i] and 'vector' in ctypes[i]:
      selector_file.write('                  '+names[i]+'new = *'+names[i]+';\n')
    elif 'pho' not in names[i] and 'vector' not in ctypes[i]:
      selector_file.write('                  '+names[i]+'new = '+names[i]+';\n')
  selector_file.write('''                  first_time = false;
                }
                phoEt_pass_count++;
              }
              phoPhi_pass_count++;
            }
            phoPFChWorstIso_pass_count++;
          }
          phoHoverE_pass_count++;
        }
        phoEta_pass_count++;
      }
      pho_total_count++;
    }
    int i_phoEt_max;
    i_phoEt_max = max_index(phoEtnew);
    if (i_phoEt_max != -1) {
''')
  for i in range(len(treelist)):# save only the photon with max Et
    if 'pho' in names[i] and 'vector' in ctypes[i]:
      selector_file.write('      '+names[i]+'new = {'+names[i]+'new.at(i_phoEt_max)};\n')
  selector_file.write('      newtree.Fill();\n')# fill the new TTree
  selector_file.write('    }\n')
  for i in range(len(treelist)):
    if 'vector' in ctypes[i]:# clear all new vectors
      selector_file.write('    '+names[i]+'new.clear();\n')
  selector_file.write('''  }
  newtree.Write();
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
  Float_t n_events_hovere_postcut;
  Float_t n_events_pfchworstiso_postcut;
  Float_t n_events_phi_postcut;
  newstatstree.Branch("n_events_precut",&n_events_precutnew,"n_events_precut/F");
  newstatstree.Branch("n_events_t_postcut",&n_events_t_postcutnew,"n_events_t_postcut/F");
  newstatstree.Branch("n_events_tm_postcut",&n_events_tm_postcutnew,"n_events_tm_postcut/F");
  newstatstree.Branch("n_events_eta_postcut",&n_events_eta_postcut,"n_events_eta_postcut/F");
  newstatstree.Branch("n_events_hovere_postcut",&n_events_hovere_postcut,"n_events_hovere_postcut/F");
  newstatstree.Branch("n_events_pfchworstiso_postcut",&n_events_pfchworstiso_postcut,"n_events_pfchworstiso_postcut/F");
  newstatstree.Branch("n_events_phi_postcut",&n_events_phi_postcut,"n_events_phi_postcut/F");
  n_events_precutnew = n_events_precut;
  n_events_t_postcutnew = n_events_t_postcut;
  n_events_tm_postcutnew = n_events_tm_postcut;
  n_events_eta_postcut = phoEta_pass_count;
  n_events_hovere_postcut = phoHoverE_pass_count;
  n_events_pfchworstiso_postcut = phoPFChWorstIso_pass_count;
  n_events_phi_postcut = phoPhi_pass_count;
  cout << "total entries: " << n << '\\n';
  cout << "prompt pass count: " << phoPhi_pass_count << '\\n';
  cout << "pass rate: " << phoPhi_pass_count/fn << '\\n';
  newstatstree.Fill();
  newstatstree.Write();
  fnew.Close();
  f.Close();
  return 0;
}''')
  selector_file.close()

def make_beamhalo_filter(dirstring,treename):
  selector_file = open("beamhalo_filter.cpp",'w')
  selector_file.write('''#include <cstdio>
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

main(int argc, char** argv) {
  string inputfilename = argv[1];
  ifstream inputfile;
  string inputcontent;
  inputfile.open(inputfilename);
  getline(inputfile,inputcontent);
  inputfile.close();
  string infile = inputcontent;
  TFile f(infile.c_str());
  TTree* tree = (TTree*)f.Get("EventTree");
  tree->SetBranchStatus("*",1);
  string outname = barename(infile)+"_promptcut.root";
  TFile fnew(outname.c_str(),"RECREATE");
  TTree newtree("EventTree","Event data");
''')
  treelist = parse_ttree(dirstring,treename)
  ctypes = []
  names = []
  for i in range(len(treelist)):
    entry = treelist[i]
    ctypes.append(entry['cpptype'])
    names.append(entry['name'])
  for i in range(len(treelist)):# infile branch initialization
    if 'vector' in ctypes[i]:
      selector_file.write('  '+ctypes[i]+'* '+names[i]+' = nullptr;\n')
    else:
      selector_file.write('  '+ctypes[i]+' '+names[i]+';\n')
  for i in range(len(treelist)):# infile SetBranchAddress
    selector_file.write('  tree->SetBranchAddress("'+names[i]+'",&'+names[i]+');\n')
  for i in range(len(treelist)):# outfile branch initialization
    selector_file.write('  '+ctypes[i]+' '+names[i]+'new;\n')
  for i in range(len(treelist)):# outfile branch initialization
    selector_file.write('  newtree.Branch("'+names[i]+'",&'+names[i]+'new);\n')
  selector_file.write('''  float prompt_pass_count = 0.0f;
  float pho_total_count = 0.0f;
  int n = tree->GetEntries();
  float fn = n;
  float phoEta_pass_count = 0.0f;
  float phoPhi_pass_count = 0.0f;
  bool first_time;
  for (Int_t i=0; i<n; i++) {
    tree->GetEntry(i);
    first_time = true;
    for (int j=0; j<phoEta->size(); j++) {
      if(fabs(phoEta->at(j)) > 1.566 && fabs(phoEta->at(j)) < 2.5) {
        if (fabs(phoPhi->at(j)) > 2.8 || fabs(phoPhi->at(j)) < 0.2 ) {
''')
  for i in range(len(treelist)):
    if 'pho' in names[i] and 'vector' in ctypes[i]:# keep photon info only for passing photons
      selector_file.write('          '+names[i]+'new.push_back('+names[i]+'->at(j));\n')
  selector_file.write('          if (first_time) {\n')
  for i in range(len(treelist)):# the first time a passing photon is found in an event, save all non-pho vectors
    if 'pho' not in names[i] and 'vector' in ctypes[i]:
      selector_file.write('            '+names[i]+'new = *'+names[i]+';\n')
    elif 'pho' not in names[i] and 'vector' not in ctypes[i]:
      selector_file.write('            '+names[i]+'new = '+names[i]+';\n')
  selector_file.write('''            first_time = false;
          }
          phoPhi_pass_count++;
        }
        phoEta_pass_count++;
      }
      pho_total_count++;
    }
    if (phoEnew.size() == 0) {
''')
  for i in range(len(treelist)):
    if 'pho' in names[i] and 'vector' in ctypes[i]:
      selector_file.write('      '+names[i]+'new.clear();\n')
  selector_file.write('''    }
    else {
      newtree.Fill();
    }
''')
  for i in range(len(treelist)):
    if 'vector' in ctypes[i]:# clear all new vectors
      selector_file.write('    '+names[i]+'new.clear();\n')
  selector_file.write('''  }
  newtree.Write();
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
  Float_t n_events_phi_postcut;
  newstatstree.Branch("n_events_precut",&n_events_precutnew,"n_events_precut/F");
  newstatstree.Branch("n_events_t_postcut",&n_events_t_postcutnew,"n_events_t_postcut/F");
  newstatstree.Branch("n_events_tm_postcut",&n_events_tm_postcutnew,"n_events_tm_postcut/F");
  newstatstree.Branch("n_events_eta_postcut",&n_events_eta_postcut,"n_events_eta_postcut/F");
  newstatstree.Branch("n_events_phi_postcut",&n_events_phi_postcut,"n_events_phi_postcut/F");
  n_events_precutnew = n_events_precut;
  n_events_t_postcutnew = n_events_t_postcut;
  n_events_tm_postcutnew = n_events_tm_postcut;
  n_events_eta_postcut = phoEta_pass_count;
  n_events_phi_postcut = phoPhi_pass_count;
  cout << "total entries: " << n << '\\n';
  cout << "prompt pass count: " << phoPhi_pass_count << '\\n';
  cout << "pass rate: " << phoPhi_pass_count/fn << '\\n';
  newstatstree.Fill();
  newstatstree.Write();
  fnew.Close();
  f.Close();
  return 0;
}''')

phoR9_d = {
    "name":"phoR9",
    "low":0,
    "high":1,
    "units":"R9",
    "position":0}
phoSigmaIEtaIEtaFull5x5_d = {
    "name":"phoSigmaIEtaIEtaFull5x5",
    "low":0,
    "high":0.07,
    "units":"\sigma_{i\eta i\eta}",
    "position":1}
phoSigmaIEtaIPhiFull5x5_d = {
    "name":"phoSigmaIEtaIPhiFull5x5",
    "low":-0.002,
    "high":0.002,
    "units":"\sigma_{i\eta i\phi}",
    "position":2}
phoSigmaIPhiIPhiFull5x5_d = {
    "name":"phoSigmaIPhiIPhiFull5x5",
    "low":0,
    "high":0.1,
    "units":"\sigma_{i\phi i\phi}",
    "position":3}
phoSeedTime_d = {
    "name":"phoSeedTime",
    "low":-140,
    "high":140,
    "units":"SeedTime",
    "position":4}
#phoPhi_d = {
#    "name":"phoPhi",
#    "low":-3.1416,
#    "high":3.1416,
#    "units":"\phi",
#    "position":5}
phoHoverE_d = {
    "name":"phoHoverE",
    "low":0,
    "high":5,
    "units":"H/E",
    "position":5}
phoNumHERH_d = {
    "name":"phoNumHERH",
    "low":0,
    "high":10,
    "units":"num HE rec-hits",
    "position":6}
phoSCEtaWidth_d = {
    "name":"phoSCEtaWidth",
    "low":0,
    "high":5,
    "units":"pho SC Eta Width",
    "position":7}
phoSCPhiWidth_d = {
    "name":"phoSCPhiWidth",
    "low":0,
    "high":5,
    "units":"pho SC Phi Width",
    "position":8}
var_od = OrderedDict()
var_od['phoR9'] = phoR9_d
var_od['phoSigmaIEtaIEtaFull5x5'] = phoSigmaIEtaIEtaFull5x5_d
var_od['phoSigmaIEtaIPhiFull5x5'] = phoSigmaIEtaIPhiFull5x5_d
var_od['phoSigmaIPhiIPhiFull5x5'] = phoSigmaIPhiIPhiFull5x5_d
var_od['phoSeedTime'] = phoSeedTime_d
#var_od['phoPhi'] = phoPhi_d
var_od['phoHoverE'] = phoHoverE_d
var_od['phoNumHERH'] = phoNumHERH_d
var_od['phoSCEtaWidth'] = phoSCEtaWidth_d
var_od['phoSCPhiWidth'] = phoSCPhiWidth_d

def make_corr(dirstring,treename):
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_treepath(dirstring,treename))
  r.gStyle.SetOptStat(0)
  #r.gStyle.SetTitleAlign(33)
  #r.gStyle.SetTitleX(0.5)
  hcorr = r.TH2F('hcorr','Pearson Correlation Matrix '+get_barename(get_filestring(dirstring)),len(var_od),0,len(var_od),len(var_od),0,len(var_od))
  for var1 in var_od:
    varx = var_od[var1]
    hcorr.GetXaxis().SetBinLabel(varx['position']+1,varx['units'])
    hcorr.GetYaxis().SetBinLabel(varx['position']+1,varx['units'])
    for var2 in var_od:
      vary = var_od[var2]
      h = r.TH2F('h',varx['name']+" vs "+vary['name'],100,varx['low'],varx['high'],100,vary['low'],vary['high'])
      h.GetXaxis().SetTitle(varx['units'])
      h.GetYaxis().SetTitle(vary['units'])
      #sel = raw_input("Selection (continue for none): ")
      t.Draw(vary['name']+":"+varx['name']+">>h",'','goff')
      corfac = h.GetCorrelationFactor()
      hcorr.Fill(varx['position'],vary['position'],corfac)
      #print "Pearson Correlation Factor: "+str(corfac)
  canv = r.TCanvas('corr_canvas','Pearson Correlation Matrix')
  hcorr.Draw('col text')
  #canv.Update()
  save = raw_input("Save? (y/n): ")
  if save == 'y':
    canv.SaveAs('corr_matrix_'+get_barename(get_filestring(dirstring))+'.png')

def process_ttree(dirstring,treename):
  print "TTree: "+treename
  chosen = {}
  print "c: compare branch from 2 files"
  print "1: histogram a single branch"
  print "2: histogram 2 branches"
  #print "pf: make a prompt filter"
  #print "bhf: make a beam halo filter"
  print "corr: make a correlation matrix"
  action = raw_input("Please choose an action: ")
  if action == '1':
    plot_treebranch(dirstring,treename)
  elif action == 'c':
    compare_var(dirstring,treename)
  elif action == '2':
    process_x_v_y(dirstring,treename)
  elif action == 'sqrt':
    plot_he_r(dirstring,treename)
  elif action == 'phoSC':
    plot_HE_phoSC_prox(dirstring,treename)
  elif action == 'pf':
    make_prompt_filter(dirstring,treename)
  #elif action == 'bhf':
  #  make_beamhalo_filter(dirstring,treename)
  elif action == 'corr':
    make_corr(dirstring,treename)
  return chosen

def process_TH1(dirstring,th1name):
  print "Got "+th1string
  # What do I want to do with a histogram?
  # I want to display it.

def process(c):
  chosen = {}
  dirstring = c['dir']
  ctype = c['cpptype']
  name = c['name']
  if ctype == 'TDirectoryFile':
    chosen = process_directory(dirstring+':'+c['name'])
  elif ctype == 'TTree':
    process_ttree(dirstring,c['name'])
  elif 'TH1' in ctype:
    process_TH1(dirstring,c['name'])
  return chosen

c0 = process_directory(filestring)
c1 = process(c0)
while c1 != {}:
  c1 = process(c1)
