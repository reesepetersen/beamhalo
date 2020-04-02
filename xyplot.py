#!/usr/bin/env python
import os
import sys
import numpy as np
from ROOT import *
import datetime
from array import array
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.colors import ListedColormap
import matplotlib as mpl

def get_barename(filestring):
  if '/' in filestring:
    return filestring.split('/')[-1].split('.')[0]
  else:
    return filestring.split('.')[0]

def get_phiewa(rhphi,rhside,rhe):
  if len(rhe) > 0:
    return np.average(rhphi[rhside],weights=rhe)
  else:
    return rhe

def get_pointer_xy(rhphi,rhside,rhe):
  if len(rhe) > 0:
    prad = 30
    rhphiav = np.average(rhphi[rhside],weights=rhe)
    rhx = np.linspace(0,prad*np.cos(rhphiav),20)
    rhy = np.linspace(0,prad*np.sin(rhphiav),20)
    return rhx, rhy
  else:
    return rhe, rhe

infile = sys.argv[1]
rf = TFile.Open(infile)
t = rf.Get("EventTree")
Nentries = t.GetEntries()
ps_rmax = 127#cm
he_rmax = 169#cm
shower_max_z = 319.5#cm
do_pho = False
auto_max = True
en_weighted_phi = True
if do_pho:
  en_weighted_phi = False
#for iEntry in range(Nentries):
for iEntry in range(25):
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
  if do_pho:
    phosceta = np.array(getattr(t,'phoSCEta'))
    phoscphi = np.array(getattr(t,'phoSCPhi'))
    phosce = np.array(getattr(t,'phoSCE'))
    phoscr = shower_max_z/np.sinh(phosceta)
    phoscx = phoscr*np.cos(phoscphi)
    phoscy = phoscr*np.sin(phoscphi)
    phoscpos = phosceta > 1.65
    phoscneg = phosceta < -1.65
    phoscxpos = phoscx[phoscpos]
    phoscxneg = phoscx[phoscneg]
    phoscypos = phoscy[phoscpos]
    phoscyneg = phoscy[phoscneg]
    phoscepos = phosce[phoscpos]
    phosceneg = phosce[phoscneg]
    bpos = plt.hist2d(phoscxpos, phoscypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens, weights=phoscepos)
    bneg = plt.hist2d(phoscxneg, phoscyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens, weights=phosceneg)
  else:
    eerheta = np.array(getattr(t,'eerhEta'))
    eerhe = np.array(getattr(t,'eerhE'))
    eerhx = np.array(getattr(t,'eerhX'))
    eerhy = np.array(getattr(t,'eerhY'))
    eerhpos = eerheta > 1.65
    eerhneg = eerheta < -1.65
    eerhxpos = eerhx[eerhpos]
    eerhxneg = eerhx[eerhneg]
    eerhypos = eerhy[eerhpos]
    eerhyneg = eerhy[eerhneg]
    eerhepos = eerhe[eerhpos]
    eerheneg = eerhe[eerhneg]
    bpos = plt.hist2d(eerhxpos, eerhypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens, weights=eerhepos)
    bneg = plt.hist2d(eerhxneg, eerhyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens, weights=eerheneg)
  fig1 = plt.figure()
  if en_weighted_phi:
    esrhphi = np.array(getattr(t,'esrhPhi'))
    eerhphi = np.array(getattr(t,'eerhPhi'))
    hbherhphi = np.array(getattr(t,'hbherhPhi'))
    #esrhphiavpos = get_phiewa(esrhphi,esrhpos,esrhepos)
    #esrhphiavneg = get_phiewa(esrhphi,esrhneg,esrheneg)
    #eerhphiavpos = get_phiewa(eerhphi,eerhpos,eerhepos)
    #eerhphiavneg = get_phiewa(eerhphi,eerhneg,eerheneg)
    #hbherhphiavpos = get_phiewa(hbherhphi,hbherhpos,hbherhepos)
    #hbherhphiavneg = get_phiewa(hbherhphi,hbherhneg,hbherheneg)
    prad = 30 # p for pointer (in energy-weighted average phi)
    #esrhpxpos = np.linspace(0,prad*np.cos(esrhphiavpos),100)
    #esrhpxneg = np.linspace(0,prad*np.cos(esrhphiavneg),100)
    #esrhpypos = np.linspace(0,prad*np.sin(esrhphiavpos),100)
    #esrhpyneg = np.linspace(0,prad*np.sin(esrhphiavneg),100)
    #eerhpxpos = np.linspace(0,prad*np.cos(eerhphiavpos),100)
    #eerhpxneg = np.linspace(0,prad*np.cos(eerhphiavneg),100)
    #eerhpypos = np.linspace(0,prad*np.sin(eerhphiavpos),100)
    #eerhpyneg = np.linspace(0,prad*np.sin(eerhphiavneg),100)
    #hbherhpxpos = np.linspace(0,prad*np.cos(hbherhphiavpos),100)
    #hbherhpxneg = np.linspace(0,prad*np.cos(hbherhphiavneg),100)
    #hbherhpypos = np.linspace(0,prad*np.sin(hbherhphiavpos),100)
    #hbherhpyneg = np.linspace(0,prad*np.sin(hbherhphiavneg),100)
    esrhpxpos, esrhpypos = get_pointer_xy(esrhphi,esrhpos,esrhepos)
    esrhpxneg, esrhpyneg = get_pointer_xy(esrhphi,esrhneg,esrheneg)
    eerhpxpos, eerhpypos = get_pointer_xy(eerhphi,eerhpos,eerhepos)
    eerhpxneg, eerhpyneg = get_pointer_xy(eerhphi,eerhneg,eerheneg)
    hbherhpxpos, hbherhpypos = get_pointer_xy(hbherhphi,hbherhpos,hbherhepos)
    hbherhpxneg, hbherhpyneg = get_pointer_xy(hbherhphi,hbherhneg,hbherheneg)
    rphipos = plt.hist2d(esrhpxpos, esrhpypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Reds)
    rphineg = plt.hist2d(esrhpxneg, esrhpyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Reds)
    gphipos = plt.hist2d(eerhpxpos, eerhpypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens)
    gphineg = plt.hist2d(eerhpxneg, eerhpyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Greens)
    bphipos = plt.hist2d(hbherhpxpos, hbherhpypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Blues)
    bphineg = plt.hist2d(hbherhpxneg, hbherhpyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Blues)
    rphip = np.array(rphipos[0])
    rphin = np.array(rphineg[0])
    gphip = np.array(gphipos[0])
    gphin = np.array(gphineg[0])
    bphip = np.array(bphipos[0])
    bphin = np.array(bphineg[0])
  rpos = plt.hist2d(esrhxpos, esrhypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Reds, weights=esrhepos)
  rneg = plt.hist2d(esrhxneg, esrhyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Reds, weights=esrheneg)
  gpos = plt.hist2d(hbherhxpos, hbherhypos, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Blues, weights=hbherhepos)
  gneg = plt.hist2d(hbherhxneg, hbherhyneg, 100, ((-he_rmax,he_rmax),(-he_rmax,he_rmax)), cmap=plt.cm.Blues, weights=hbherheneg)
  rp = np.array(rpos[0])
  bp = np.array(gpos[0])
  gp = np.array(bpos[0])
  rn = np.array(rneg[0])
  bn = np.array(gneg[0])
  gn = np.array(bneg[0])
  rpmax = np.amax(rp)
  bpmax = np.amax(bp)
  gpmax = np.amax(gp)
  rnmax = np.amax(rn)
  bnmax = np.amax(bn)
  gnmax = np.amax(gn)
  ch = np.linspace(0,1,256)
  full = np.ones_like(ch)
  red = np.zeros((256,4))
  green = np.zeros((256,4))
  blue = np.zeros((256,4))
  red[:,0] = ch
  green[:,1] = ch
  blue[:,2] = ch
  red[:,3] = full
  green[:,3] = full
  blue[:,3] = full
  cmr = ListedColormap(red)
  cmg = ListedColormap(green)
  cmb = ListedColormap(blue)
  if rpmax > rnmax:
    imr = plt.pcolormesh(rp, cmap=cmr)
  elif rnmax > rpmax:
    imr = plt.pcolormesh(rn, cmap=cmr)
  if gpmax > gnmax:
    img = plt.pcolormesh(gp, cmap=cmg)
  elif gnmax > gpmax:
    img = plt.pcolormesh(gn, cmap=cmg)
  if bpmax > bnmax:
    imb = plt.pcolormesh(bp, cmap=cmb)
  elif bnmax > bpmax:
    imb = plt.pcolormesh(bn, cmap=cmb)
  if auto_max:
    rmax = max(rpmax,rnmax)
    bmax = max(bpmax,bnmax)
    gmax = max(gpmax,gnmax)
  else:
    rmax = 0.06
    bmax = 350
    gmax = 1000
  if rmax > 0.:
    rp = rp/rmax
    rn = rn/rmax
  if bmax > 0.:
    bp = bp/bmax
    bn = bn/bmax
  if gmax > 0.:
    gp = gp/gmax
    gn = gn/gmax
  rgbpos = np.dstack((rp,gp,bp))
  rgbneg = np.dstack((rn,gn,bn))
  if en_weighted_phi:
    rgbphip = np.dstack((rphip,gphip,bphip))
    rgbphin = np.dstack((rphin,gphin,bphin))
    rgbpos += rgbphip
    rgbneg += rgbphin
  pos = ndimage.rotate(rgbpos,90)
  neg = ndimage.rotate(rgbneg,90)
  plt.close(fig1)
  fig2 = plt.figure(figsize=(10,10))
  grid = plt.GridSpec(4, 6, wspace=0.3, hspace=0.2)
  box1 = plt.subplot(grid[0:3,0:3])
  box2 = plt.subplot(grid[0:3,3:6])
  event = np.array(getattr(t,'event'))
  fig2.suptitle(get_barename(infile)+':'+str(event)+' : '+str(iEntry))
  i1 = box1.imshow(pos, extent=[-he_rmax,he_rmax,-he_rmax,he_rmax], interpolation='none')
  box1.set_aspect('equal')
  box1.set_ylabel('y (cm)')
  box1.set_xlabel('x (cm)')
  box1.set_title('+z')
  i2 = box2.imshow(neg, extent=[-he_rmax,he_rmax,-he_rmax,he_rmax], interpolation='none')
  box2.set_aspect('equal')
  box2.set_xlabel('x (cm)')
  box2.set_title('-z')
  box2.set_yticklabels([])
  caxr = fig2.add_axes([0.1,0.3,0.8,0.03])
  caxg = fig2.add_axes([0.1,0.2,0.8,0.03])
  caxb = fig2.add_axes([0.1,0.1,0.8,0.03])
  normr = mpl.colors.Normalize(vmin=0, vmax=rmax)
  normg = mpl.colors.Normalize(vmin=0, vmax=gmax)
  normb = mpl.colors.Normalize(vmin=0, vmax=bmax)
  cobar = mpl.colorbar.ColorbarBase(caxr, cmap=cmr, orientation='horizontal', norm=normr)
  cobag = mpl.colorbar.ColorbarBase(caxg, cmap=cmg, orientation='horizontal', norm=normg)
  cobab = mpl.colorbar.ColorbarBase(caxb, cmap=cmb, orientation='horizontal', norm=normb)
  cobar.set_label('Energy deposited in Preshower (GeV)')
  cobag.set_label('Energy deposited in ECal End-cap (GeV)')
  cobab.set_label('Energy deposited in HCAL End-cap (GeV)')
  #plt.show()
  plt.savefig(get_barename(infile)+'_end_cap_evt_'+str(iEntry)+'.png')
