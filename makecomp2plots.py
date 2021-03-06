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
import numpy as np

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

phoR9_d = {
    "name":"phoR9",
    "low":0.4,
    "high":1,
    "units":"R9",
    "position":0,
    "yscale":"linear",
    "nbins":100}
phoSigmaIEtaIEtaFull5x5_d = {
    "name":"phoSigmaIEtaIEtaFull5x5",
    "low":0,
    "high":0.05,
    "units":"\sigma_{i\eta i\eta}",
    "position":1,
    "yscale":"linear",
    "nbins":100}
phoSigmaIEtaIPhiFull5x5_d = {
    "name":"phoSigmaIEtaIPhiFull5x5",
    "low":-0.0006,
    "high":0.0006,
    "units":"\sigma_{i\eta i\phi}",
    "position":2,
    "yscale":"linear",
    "nbins":100}
phoSigmaIPhiIPhiFull5x5_d = {
    "name":"phoSigmaIPhiIPhiFull5x5",
    "low":0,
    "high":0.07,
    "units":"\sigma_{i\phi i\phi}",
    "position":3,
    "yscale":"linear",
    "nbins":100}
phoSeedTime_d = {
    "name":"phoSeedTime",
    "low":-3.1,
    "high":3.1,
    "units":"SeedTime",
    "position":4,
    "yscale":"log",
    "nbins":100}
phoPhi_d = {
    "name":"phoPhi",
    "low":-3.1416,
    "high":3.1416,
    "units":"\phi",
    "position":5,
    "yscale":"linear",
    "nbins":100}
phoHoverE_d = {
    "name":"phoHoverE",
    "low":0,
    "high":1.4,
    "units":"H/E",
    "position":6,
    "yscale":"log",
    "nbins":100}
phoNumHERHzside_d = {
    "name":"phoNumHERHzside",
    "low":0,
    "high":12,
    "units":"num HE rec-hits",
    "position":7,
    "yscale":"linear",
    "nbins":12}
phoNumESRHzside_d = {
    "name":"phoNumESRHzside",
    "low":0,
    "high":12,
    "units":"num ES rec-hits",
    "position":8,
    "yscale":"linear",
    "nbins":12}
phoSCEtaWidth_d = {
    "name":"phoSCEtaWidth",
    "low":0,
    "high":0.03,
    "units":"pho SC Eta Width",
    "position":9,
    "yscale":"linear",
    "nbins":100}
phoSCPhiWidth_d = {
    "name":"phoSCPhiWidth",
    "low":0,
    "high":0.14,
    "units":"pho SC Phi Width",
    "position":10,
    "yscale":"linear",
    "nbins":100}
phoE_d = {
    "name":"phoE",
    "low":0,
    "high":5000,
    "units":"phoE",
    "position":11,
    "yscale":"linear",
    "nbins":100}
phoEt_d = {
    "name":"phoEt",
    "low":0,
    "high":2500,
    "units":"phoEt",
    "position":12,
    "yscale":"log",
    "nbins":100}
phoHaloHE_d = {
    "name":"phoHaloHE",
    "low":0,
    "high":100,
    "units":"Halo HE",
    "position":13,
    "yscale":"log",
    "nbins":100}
phoHaloPre_d = {
    "name":"phoHaloPre",
    "low":0,
    "high":0.45,
    "units":"Halo Pre",
    "position":14,
    "yscale":"log",
    "nbins":100}
phoEnWeAn_d = {
    "name":"phoEnWeAn",
    "low":-1.2,
    "high":1.2,
    "units":"Energy-weighted Angle",
    "position":15,
    "yscale":"log",
    "nbins":100}
phoA1_d = {
    "name":"phoA1",
    "low":-0.9,
    "high":1.2,
    "units":"Region 1 Angle",
    "position":16,
    "yscale":"log",
    "nbins":100}
phoA2_d = {
    "name":"phoA2",
    "low":-0.9,
    "high":1.2,
    "units":"Region 2 Angle",
    "position":17,
    "yscale":"log",
    "nbins":100}
phoE1_d = {
    "name":"phoE1",
    "low":0,
    "high":120,
    "units":"Region 1 Energy",
    "position":18,
    "yscale":"log",
    "nbins":100}
phoE2_d = {
    "name":"phoE2",
    "low":0,
    "high":120,
    "units":"Region 2 Energy",
    "position":19,
    "yscale":"log",
    "nbins":100}
phoPreEnWeAn_d = {
    "name":"phoPreEnWeAn",
    "low":-1,
    "high":1.5,
    "units":"Energy-weighted Angle",
    "position":15,
    "yscale":"log",
    "nbins":100}
phoPreA1_d = {
    "name":"phoPreA1",
    "low":-1,
    "high":1.5,
    "units":"Region 1 Angle",
    "position":16,
    "yscale":"log",
    "nbins":100}
phoPreA2_d = {
    "name":"phoPreA2",
    "low":-1,
    "high":1.5,
    "units":"Region 2 Angle",
    "position":17,
    "yscale":"log",
    "nbins":100}
phoPreE1_d = {
    "name":"phoPreE1",
    "low":0,
    "high":3,
    "units":"Region 1 Energy",
    "position":18,
    "yscale":"log",
    "nbins":100}
phoPreE2_d = {
    "name":"phoPreE2",
    "low":0,
    "high":3,
    "units":"Region 2 Energy",
    "position":19,
    "yscale":"log",
    "nbins":100}
var_od = OrderedDict()
#var_od['phoR9'] = phoR9_d
#var_od['phoSigmaIEtaIEtaFull5x5'] = phoSigmaIEtaIEtaFull5x5_d
#var_od['phoSigmaIEtaIPhiFull5x5'] = phoSigmaIEtaIPhiFull5x5_d
#var_od['phoSigmaIPhiIPhiFull5x5'] = phoSigmaIPhiIPhiFull5x5_d
#var_od['phoSeedTime'] = phoSeedTime_d
#var_od['phoPhi'] = phoPhi_d
#var_od['phoHoverE'] = phoHoverE_d
#var_od['phoNumHERHzside'] = phoNumHERHzside_d
#var_od['phoNumESRHzside'] = phoNumESRHzside_d
#var_od['phoSCEtaWidth'] = phoSCEtaWidth_d
#var_od['phoSCPhiWidth'] = phoSCPhiWidth_d
#var_od['phoE'] = phoE_d
#var_od['phoEt'] = phoEt_d
#var_od['phoHaloHE'] = phoHaloHE_d
#var_od['phoHaloPre'] = phoHaloPre_d
#var_od['phoEnWeAn'] = phoEnWeAn_d
#var_od['phoE1'] = phoE1_d
#var_od['phoE2'] = phoE2_d
#var_od['phoA1'] = phoA1_d
#var_od['phoA2'] = phoA2_d
var_od['phoPreEnWeAn'] = phoPreEnWeAn_d
var_od['phoPreE1'] = phoPreE1_d
var_od['phoPreE2'] = phoPreE2_d
var_od['phoPreA1'] = phoPreA1_d
var_od['phoPreA2'] = phoPreA2_d

# check input file for .root extension
infilelist = sys.argv
for f in infilelist[1:]:
#filestring = sys.argv[1]
#compfilestring = sys.argv[2]
  if f.split('.')[-1] != 'root':
    print "Input file {0} does not have .root extentsion. exiting...".format(f)
    quit()

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

def get_objectpath(dirstring,obj,joiner = '/'):
  if ':' in dirstring:
    dirlist = dirstring.split(':')
    dirlist.pop(0)
    return joiner.join(dirlist)+joiner+obj
  else:
    return obj

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
    tree = rootfile.Get(get_objectpath(dirstring,treename))
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

def plot_treebranch(dirstring,treename,sel):
  choice = select_from_ttree(dirstring,treename)
  rootfile = r.TFile(get_filestring(dirstring))
  tree = rootfile.Get(get_objectpath(dirstring,treename))
  c = r.TCanvas('phoPhiMIPTotalcanvas',sel)
  if sel == "":
    sel = raw_input("Selection (default none): ")
  tree.Draw(choice['name'],sel)
  raw_input("Enter to quit ")

def plot_treebranch_MIP(dirstring,treename,sel):
  rootfile = r.TFile(get_filestring(dirstring))
  tree = rootfile.Get(get_objectpath(dirstring,treename))
  c = r.TCanvas('phoPhiMIPTotalcanvas',sel)
  tree.Draw("phoPhi",sel)
  c.SaveAs("MIPTotal/phoPhi_MIPTotal_gt_"+sel.split('>')[-1].split('&&')[0]+get_barename(get_filestring(dirstring))+".png")

def plot_treebranch_AHE(dirstring,treename,sel):
  rootfile = r.TFile(get_filestring(dirstring))
  tree = rootfile.Get(get_objectpath(dirstring,treename))
  c = r.TCanvas('phoPhiAHETotalcanvas',sel)
  tree.Draw("phoPhi",sel)
  c.SaveAs("AHETotal/phoPhi_AHETotal_gt_"+sel.split('>')[1].split('&&')[0]+get_barename(get_filestring(dirstring))+".png")

def get_color(filestring):
  if 'prompt' in filestring and 'data' in filestring:
    return r.kRed
  elif 'prompt' in filestring and 'MC' in filestring:
    return r.kOrange+1
  elif 'beamhalo' in filestring and 'data' in filestring:
    if 'not' in filestring:
      return r.kViolet+3
    else:
      return r.kBlue
  elif 'beamhalo' in filestring and 'MC' in filestring:
    return r.kViolet
  elif 'newbh' in filestring and 'MC' in filestring:
    return r.kBlue
  elif 'tmne' in filestring and 'data' in filestring:
    return r.kAzure+1
  elif 'tmn' in filestring and 'MC' in filestring:
    return r.kViolet+1
  elif 'tmn' in filestring and 'data' in filestring:
    return r.kGreen+1
  elif 'tm' in filestring and 'MC' in filestring:
    return r.kGray+2
  elif 'tm' in filestring and 'data' in filestring:
    return r.kBlue
  elif 'bhvertical' in filestring:
    return r.kCyan+1
  elif 'MC' in filestring:
    return r.kBlack

def compare_var_from_files(name,dirstring,treename,infilelist,hrlow,hrhigh,hist_units,nbins):
  file1 = infilelist[1]
  file2 = infilelist[2]
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf1 = r.TFile(file1)
  rf2 = r.TFile(file2)
  t1 = rf1.Get(get_objectpath(dirstring,treename))
  t2 = rf2.Get(get_objectpath(dirstring,treename))
  h1 = r.TH1D('h1',name,nbins,hrlow,hrhigh)
  h2 = r.TH1D('h2',name,nbins,hrlow,hrhigh)
  h1.SetLineColor(get_color(file1))
  h2.SetLineColor(get_color(file2))
  h1.GetXaxis().SetTitle(hist_units)
  sel = ''
  t1.Draw(name+'>>h1',sel,'goff')
  t2.Draw(name+'>>h2',sel,'goff')
  comp_canvas = r.TCanvas('comp_canvas','compare '+name)
  maxall = max(h1.GetMaximum(),h2.GetMaximum())
  sum1 = h1.Integral()
  sum2 = h2.Integral()
  maxsum = max(sum1,sum2)
  h1.Scale(1/sum1)
  h2.Scale(1/sum2)
  var = var_od[name]
  if var['yscale'] == "log":
    comp_canvas.SetLogy()
    h1.SetMinimum(1./(2*maxsum))
    h2.SetMinimum(1./(2*maxsum))
  max1 = h1.GetMaximum()
  max2 = h2.GetMaximum()
  maxlist = [max1,max2]
  if max1 == max(maxlist):
    h1.Draw('hist')
    h2.Draw('hist same')
  elif max2 == max(maxlist):
    h2.Draw('hist')
    h1.Draw('hist same')
  leg = r.TLegend(0.5,0.9,0.9,1.0)
  leg.AddEntry(h1,get_barename(file1),'l')
  leg.AddEntry(h2,get_barename(file2),'l')
  leg.Draw('hist')
  comp_canvas.Update()
  comp_canvas.SaveAs('../slides_beamhalo/'+name+'_'+get_barename(file1)+'_'+get_barename(file2)+'.png')

def plot_2D(x_var,x_range,x_units,y_var,y_range,y_units,dirstring,treename):
  xr = x_range.split(',')
  yr = y_range.split(',')
  xlow = float(xr[0])
  ylow = float(yr[0])
  xhigh = float(xr[1])
  yhigh = float(yr[1])
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_objectpath(dirstring,treename))
  h = r.TH2F('h',x_var+" vs "+y_var,180,xlow,xhigh,180,ylow,yhigh)
  h.GetXaxis().SetTitle(x_units)
  h.GetYaxis().SetTitle(y_units)
  sel = raw_input("Selection (continue for none): ")
  t.Draw(y_var+":"+x_var+">>h",sel,'goff')
  canv = r.TCanvas('canvas2D','2D '+x_var+' vs '+y_var)
  corfac = h.GetCorrelationFactor()
  print "Pearson Correlation Factor: "+str(corfac)
  h.SetTitle(sel)
  h.Draw('colz')
  canv.Update()
  save = raw_input("Save? (y/n): ")
  if save == 'y':
    canv.SaveAs(x_var+'_vs_'+y_var+'_'+get_barename(get_filestring(dirstring))+'.png')

def plot_2d(dirstring,treename,feat,feat_od,outdir):
  r.gStyle.SetOptStat(0)
  for ofeatname in feat_od:
    ofeat = feat_od[ofeatname]
    x = feat['name']
    y = ofeat['name']
    xnbins = feat['nbins']
    ynbins = ofeat['nbins']
    xlow = feat['low']
    ylow = ofeat['low']
    xhigh = feat['high']
    yhigh = ofeat['high']
    xunits = feat['units']
    yunits = ofeat['units']
    rf = r.TFile(get_filestring(dirstring))
    t = rf.Get(get_objectpath(dirstring,treename))
    h = r.TH2D('h2d',y+' vs '+x,xnbins,xlow,xhigh,ynbins,ylow,yhigh)
    h.GetXaxis().SetTitle(xunits)
    h.GetYaxis().SetTitle(yunits)
    sel = 'fabs(phoEta)>=1.65&&fabs(phoEta)<1.8&&phoE<5000&&phoEt>50&&pfMET>200&&HLTPho>>9&1'
    t.Draw(y+':'+x+'>>h2d',sel)
    can = r.TCanvas('can_'+y+'_vs_'+x,'2D '+x+' vs '+y)
    h.SetTitle(sel)
    h.Draw('colz')
    if 'phoHalo' in x or 'phoHalo' in y:
      can.SetLogz()
    can.Update()
    can.SaveAs(outdir+'/'+x+'_vs_'+y+'_'+get_barename(get_filestring(dirstring))+'.png')

def plot_2d_auto(dirstring,treename):
  outdir = raw_input('output dir: ')
  os.makedirs(outdir)
  feat_od = var_od
  for featname in feat_od:
    feat = feat_od[featname]
    feat_od.pop(featname)
    plot_2d(dirstring,treename,feat,feat_od,outdir)

def compare_vars(dirstring,treename,infilelist):
  for varname in var_od:
    var = var_od[varname]
    compare_var_from_files(varname,dirstring,treename,infilelist,var['low'],var['high'],var['units'],var['nbins'])

def process_x_v_y(dirstring,treename):
#  choice = select_from_ttree(dirstring,treename)
#  x_var = choice['name']
#  print "x variable: "+x_var
  x_var = raw_input("x variable: ")
  x_range = raw_input("Please give the x range (e.g. -3.14,3.14): ")
  x_units = raw_input("Please give the x units: ")
  y_var = raw_input("y variable: ")
  y_range = raw_input("Please give the y range (e.g. -3.14,3.14): ")
  y_units = raw_input("Please give the y units: ")
  plot_2D(x_var,x_range,x_units,y_var,y_range,y_units,dirstring,treename)

def make_corr(dirstring,treename):
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_objectpath(dirstring,treename))
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
  canv = r.TCanvas('corr_canvas','Pearson Correlation Matrix')
  hcorr.Draw('col text')
  #canv.Update()
  save = raw_input("Save? (y/n): ")
  if save == 'y':
    canv.SaveAs('corr_matrix_'+get_barename(get_filestring(dirstring))+'.png')

def make_reweight(dirstring,treename):
  name = "phoEt"
  var = var_od[name]
  hrlow = var['low']
  hrhigh = var['high']
  hist_units = var['units']
  nbins = var['nbins']
  file1 = infilelist[1]
  file2 = infilelist[2]
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf1 = r.TFile(file1)
  rf2 = r.TFile(file2,"UPDATE")
  t1 = rf1.Get(get_objectpath(dirstring,treename))
  t2 = rf2.Get(get_objectpath(dirstring,treename))
  h1 = r.TH1D('h1',name,nbins,hrlow,hrhigh)
  h2 = r.TH1D('h2',name,nbins,hrlow,hrhigh)
  h1.SetLineColor(get_color(file1))
  h2.SetLineColor(get_color(file2))
  h1.GetXaxis().SetTitle(hist_units)
  sel = ''
  t1.Draw(name+'>>h1',sel,'goff')
  t2.Draw(name+'>>h2',sel,'goff')
  comp_canvas = r.TCanvas('comp_canvas','compare '+name)
  hw = h1.Clone('hwpt')
  hw.Divide(h2)
  hw.Draw()
  hw.Write()

def compare_reweighted(name,dirstring,treename,infilelist,hrlow,hrhigh,hist_units,nbins):
  file1 = infilelist[1]
  file2 = infilelist[2]
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetTitleX(0.5)
  rf1 = r.TFile(file1)
  rf2 = r.TFile(file2)
  t1 = rf1.Get(get_objectpath(dirstring,treename))
  t2 = rf2.Get(get_objectpath(dirstring,treename))
  h1 = r.TH1D('h1',name,nbins,hrlow,hrhigh)
  h2 = r.TH1D('h2',name,nbins,hrlow,hrhigh)
  h1.SetLineColor(get_color(file1))
  h2.SetLineColor(get_color(file2))
  h1.GetXaxis().SetTitle(hist_units)
  sel = ''
  t1.Draw(name+'>>h1',sel,'goff')
  t2.Draw(name+'>>h2','phopTweight','goff')
  comp_canvas = r.TCanvas('comp_canvas','compare '+name)
  maxall = max(h1.GetMaximum(),h2.GetMaximum())
  sum1 = h1.Integral()
  sum2 = h2.Integral()
  maxsum = max(sum1,sum2)
  h1.Scale(1/sum1)
  h2.Scale(1/sum2)
  var = var_od[name]
  if var['yscale'] == "log":
    comp_canvas.SetLogy()
    h1.SetMinimum(1./(2*maxsum))
    h2.SetMinimum(1./(2*maxsum))
  max1 = h1.GetMaximum()
  max2 = h2.GetMaximum()
  maxlist = [max1,max2]
  if max1 == max(maxlist):
    h1.Draw('hist')
    h2.Draw('hist same')
  elif max2 == max(maxlist):
    h2.Draw('hist')
    h1.Draw('hist same')
  leg = r.TLegend(0.5,0.9,0.9,1.0)
  leg.AddEntry(h1,get_barename(file1),'l')
  leg.AddEntry(h2,get_barename(file2),'l')
  leg.Draw('hist')
  comp_canvas.Update()
  comp_canvas.SaveAs('../slides_beamhalo/'+name+'_'+get_barename(file1)+'_'+get_barename(file2)+'_pTweight.png')

def compare_vars_weighted(dirstring,treename,infilelist):
  for varname in var_od:
    var = var_od[varname]
    compare_reweighted(varname,dirstring,treename,infilelist,var['low'],var['high'],var['units'],var['nbins'])

def plot_phi_MIPTotal(dirstring,treename):
  for i in np.arange(0.9,10.9,1.0):
    plot_treebranch_MIP(dirstring,treename,"phoMIPTotEnergy>"+str(round(i,2))+"&&abs(phoEta)<1.566")

def plot_phi_AHETotal(dirstring,treename):
  for i in [-0.1,0.0,0.1,1.0,5.0,10.0,50.0,100.0]:
    plot_treebranch_AHE(dirstring,treename,"phoAHETotal>"+str(round(i,3))+"&&fabs(phoEta)>1.566")

def xy(dirstring,treename):
  r.gStyle.SetOptStat(0)
  x = 'esrhX'
  y = 'esrhY'
  xnbins = 100
  ynbins = 100
  xlow = -127
  ylow = -127
  xhigh = 127
  yhigh = 127
  xunits = 'x (cm)'
  yunits = 'y (cm)'
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_objectpath(dirstring,treename))
  h = r.TH2D('h2d',y+' vs '+x,xnbins,xlow,xhigh,ynbins,ylow,yhigh)
  h.GetXaxis().SetTitle(xunits)
  h.GetYaxis().SetTitle(yunits)
  h.SetContour(2)
  sel = 'fabs(phoEta)>=1.65&&fabs(phoEta)<1.8&&phoE<5000&&phoEt>50&&pfMET>200&&HLTPho>>9&1'
  t.Draw(y+':'+x+'>>h2d',sel)
  can = r.TCanvas('can_'+y+'_vs_'+x,'2D '+x+' vs '+y)
  h.SetTitle(sel)
  h.Draw('colz')
  #if 'phoHalo' in x or 'phoHalo' in y:
  #  can.SetLogz()
  can.Update()
  can.SaveAs(x+'_vs_'+y+'_'+get_barename(get_filestring(dirstring))+'.png')

def polar_phoPhi(dirstring,treename):
  r.gStyle.SetOptStat(0)
  r.gStyle.SetTitleAlign(33)
  r.gStyle.SetGridStyle(0)
  pi = 3.14159265358979
  nbins=100
  rf = r.TFile(get_filestring(dirstring))
  t = rf.Get(get_objectpath(dirstring,treename))
  h = r.TH1D('hphi','pho phi',nbins,-pi,pi)
  t.Draw('phoPhi>>hphi','phoEt>50&&phoE<5000&&fabs(phoEta)>1.65&&fabs(phoEta)<1.8','goff')
  rad = []
  phi = []
  delta = 2*pi/nbins
  ran = np.arange(-pi+delta/2,pi,delta)
  for i in ran:
    phi.append(i)
    rad.append(h.GetBinContent(h.GetXaxis().FindBin(i)))
  pol = r.TGraphPolar(len(rad),np.asarray(phi),np.asarray(rad))
  can = r.TCanvas('can_pho_polar_phi','Photon Phi Polar Plot',600,600)
  can.SetGrid(-1,-1)
  pol.SetTitle('Photon\ \phi\ \ ')
  pol.SetLineColor(2)
  pol.SetLineWidth(2)
  pol.GetXaxis().SetLabelSize(0.01)
  pol.Draw('L')
  can.Update()
  can.SaveAs('polarPhi_'+get_barename(get_filestring(dirstring))+'.png')
  raw_input('wait')

def process_ttree(dirstring,treename):
  print "TTree: "+treename
  chosen = {}
  print "comp: compare branch from 2 files"
  print "1: histogram a single branch"
  print "2: histogram 2 branches"
  print "corr: make a correlation matrix"
  print "mrpt: make reweight pt TH1D"
  print "compr: compare a branch between 2 files, second is reweighted"
  print "phimip: plot phi with various MIPTotal cuts"
  print "phiahe: plot phi with various AHETotal cuts"
  print "2d: make 2d histograms between all chosen features"
  print "xy: make x,y 2d histograms of esrh and herh with photons and electrons"
  print "pol: make polar phi plot of ph photons"
  action = raw_input("Please choose an action: ")
  if action == '1':
    plot_treebranch(dirstring,treename)
  elif action == 'comp':
    compare_vars(dirstring,treename,infilelist)
  elif action == '2':
    process_x_v_y(dirstring,treename)
  elif action == 'corr':
    make_corr(dirstring,treename)
  elif action == 'mrpt':
    make_reweight(dirstring,treename)
  elif action == 'compr':
    compare_vars_weighted(dirstring,treename,infilelist)
  elif action == 'phimip':
    plot_phi_MIPTotal(dirstring,treename)
  elif action == 'phiahe':
    plot_phi_AHETotal(dirstring,treename)
  elif action == '2d':
    plot_2d_auto(dirstring,treename)
  elif action == 'xy':
    xy(dirstring,treename)
  elif action == 'pol':
    polar_phoPhi(dirstring,treename)
  return chosen

def process_TH1(dirstring,th1name):
  print "Got "+th1string
  f = r.TFile(get_filestring(dirstring))
  th1 = f.Get(get_objectpath(dirstring,th1name))
  th1.Draw()

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

c0 = process_directory(infilelist[1])
c1 = process(c0)
while c1 != {}:
  c1 = process(c1)
