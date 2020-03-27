#!/usr/bin/env python
import os
import sys
import numpy as np
#import xgboost as xgb
#import pickle as pk
# from xgboost import XGBClassifier, Booster
#from sklearn.preprocessing import LabelEncoder
from ROOT import *
import datetime
from array import array
import matplotlib.pyplot as plt
from scipy import ndimage

# input file
infile = sys.argv[1]
rf = TFile.Open(infile)
t = rf.Get("EventTree")
Nentries = t.GetEntries()
ps_rmax = 127#cm
he_rmax = 169#cm
shower_max_z = 319.5#cm
#for iEntry in range(Nentries):
for iEntry in range(1):
  t.GetEntry(iEntry)
  esrhx = np.array(getattr(t,'esrhX'))
  esrhy = np.array(getattr(t,'esrhY'))
  esrheta = np.array(getattr(t,'esrhEta'))
  esrhe = np.array(getattr(t,'esrhE'))
  esrhpos = esrheta > 0
  esrhneg = esrheta < 0
  esrhxpos = esrhx[esrhpos]
  esrhxneg = esrhx[esrhneg]
  esrhypos = esrhy[esrhpos]
  esrhyneg = esrhy[esrhneg]
  esrhepos = esrhe[esrhpos]
  esrheneg = esrhe[esrhneg]
  hbherhx = np.array(getattr(t,'hbherhX'))
  hbherhy = np.array(getattr(t,'hbherhY'))
  hbherheta = np.array(getattr(t,'hbherhEta'))
  hbherhe = np.array(getattr(t,'hbherhE'))
  hbherhpos = hbherheta > 1.65
  hbherhneg = hbherheta < -1.65
  hbherhxpos = hbherhx[hbherhpos]
  hbherhxneg = hbherhx[hbherhneg]
  hbherhypos = hbherhy[hbherhpos]
  hbherhyneg = hbherhy[hbherhneg]
  hbherhepos = hbherhe[hbherhpos]
  hbherheneg = hbherhe[hbherhneg]
  phosceta = np.array(getattr(t,'phoSCEta'))
  phoscphi = np.array(getattr(t,'phoSCPhi'))
  phosce = np.array(getattr(t,'phoSCE'))
  phoscr = shower_max_z/np.sinh(phosceta)
  phoscx = phoscr*np.cos(phoscphi)
  phoscy = phoscr*np.sin(phoscphi)
  phoscpos = phosceta > 1.65
  phoscneg = phosceta < 1.65
  phoscxpos = phoscx[phoscpos]
  phoscxneg = phoscx[phoscneg]
  phoscypos = phoscy[phoscpos]
  phoscyneg = phoscy[phoscneg]
  phoscepos = phosce[phoscpos]
  phosceneg = phosce[phoscneg]
  fig1 = plt.figure(1)
  espos = plt.hist2d(esrhxpos, esrhypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Reds, weights=esrhepos)
  esneg = plt.hist2d(esrhxneg, esrhyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Reds, weights=esrheneg)
  hepos = plt.hist2d(hbherhxpos, hbherhypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Blues, weights=hbherhepos)
  heneg = plt.hist2d(hbherhxneg, hbherhyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Blues, weights=hbherheneg)
  pscpos = plt.hist2d(phoscxpos, phoscypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens, weights=phoscepos)
  pscneg = plt.hist2d(phoscxneg, phoscyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens, weights=phosceneg)
  rp = np.array(espos[0])
  bp = np.array(hepos[0])
  gp = np.array(pscpos[0])
  rpmax = np.amax(rp)
  bpmax = np.amax(bp)
  gpmax = np.amax(gp)
  if rpmax > 0.:
    rp = rp/rpmax
  if bpmax > 0.:
    bp = bp/bpmax
  if gpmax > 0.:
    gp = gp/gpmax
  rn = np.array(esneg[0])
  bn = np.array(heneg[0])
  gn = np.array(pscneg[0])
  rnmax = np.amax(rn)
  bnmax = np.amax(bn)
  gnmax = np.amax(gn)
  if rnmax > 0.:
    rn = rn/rnmax
  if bnmax > 0.:
    bn = bn/bnmax
  if gnmax > 0.:
    gn = gn/gnmax
  rgbpos = np.dstack((rp,gp,bp))
  rgbneg = np.dstack((rn,gn,bn))
  pos = ndimage.rotate(rgbpos,90)
  neg = ndimage.rotate(rgbneg,90)
  plt.close(fig1)
  fig2, (box1, box2) = plt.subplots(1, 2, gridspec_kw={'wspace':0.1})
  fig2.suptitle('Preshower (red) and HE (blue) rec-hits, +z left, -z right')
  box1.imshow(pos, extent=[-he_rmax,he_rmax,-he_rmax,he_rmax], interpolation='none')
  box1.set_aspect('equal')
  box1.set_ylabel('y (cm)')
  box1.set_xlabel('x (cm)')
  box2.imshow(neg, extent=[-he_rmax,he_rmax,-he_rmax,he_rmax], interpolation='none')
  box2.set_aspect('equal')
  box2.set_xlabel('x (cm)')
  plt.show()
  #plt.savefig('end_cap_evt_'+str(iEntry)+'.png')
