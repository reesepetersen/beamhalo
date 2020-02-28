#!/usr/bin/env python
import os
import sys
import numpy as np
import xgboost as xgb
import pickle as pk
# from xgboost import XGBClassifier, Booster
from sklearn.preprocessing import LabelEncoder
from ROOT import *
import datetime
from array import array

base_directory = '/local/cms/user/pet00831/beamhalo/train_results'
output_directory = base_directory + '/result_trees/'
model_file = base_directory + '/antgc_beamhalo_bdt.pickle'
sample_list = sys.argv[1]

feature_names = ['phoR9','phoSigmaIEtaIEtaFull5x5','phoSigmaIEtaIPhiFull5x5','phoSigmaIPhiIPhiFull5x5','phoHoverE','phoSCEtaWidth','phoSCPhiWidth']
#feature_names = ['phoR9','phoSigmaIEtaIEtaFull5x5','phoSigmaIEtaIPhiFull5x5','phoSigmaIPhiIPhiFull5x5','phoHoverE','phoSeedTime','phoNumHERHzside','phoNumESRHzside','phoSCEtaWidth','phoSCPhiWidth']

def runBDTonTree(inputTreePath, XGBClf, outTreeDir):
	print 'Running BDT on tree ' + inputTreePath
	# load input tree
	inputFile = TFile.Open(inputTreePath, "read")
        inputTree = inputFile.Get("EventTree")

	outFilePath = outTreeDir + os.path.basename(inputTreePath)
	outFilePath = outFilePath.replace("EventTree", "BDT")

	# create result tree
	outFile = TFile(outFilePath, 'recreate')
	#outTree = TTree( 'BDTtree', 'BDT results' )
	outTree = inputTree.CloneTree(0);
	outTree.SetDirectory(outFile.GetDirectory(''))

	# create new branches for result tree
        _iEntry = array( 'i', [0])  
        bdtResult = array( 'i', [0])
        sigProb = array( 'f', [0.])
        bgProb = array( 'f', [0.])

        outTree.Branch("sigProb", sigProb, "sigProb/F")
        outTree.Branch("bgProb", bgProb, "bgProb/F")
        outTree.Branch("bdtResult", bdtResult, "bdtResult/I")
        outTree.Branch("iEntry", _iEntry, "iEntry/L") 
        
        outTree.SetBranchStatus("sigProb",1)
        outTree.SetBranchStatus("bgProb",1)
        outTree.SetBranchStatus("bdtResult",1)
        outTree.SetBranchStatus("iEntry",1)

	Nentries = inputTree.GetEntries()
	fracNentries = Nentries/20;

	print 'Entries in tree : %s'%(Nentries)

        XGBClf.feature_names = feature_names

	for iEntry in range(Nentries):
	#for iEntry in range(1,5000):
		if( iEntry % fracNentries == 0): print "Analyzing entry %s"%(iEntry)
		inputTree.GetEntry(iEntry)
		
		xFeats = []

		for xFeat in feature_names: xFeats.append(getattr(inputTree, xFeat))

		xFeats = np.array(xFeats, dtype=float)
                #print xFeats.shape
                xFeats = np.transpose(xFeats)
		#xFeats = np.atleast_2d(xFeats)
                
                _bdtResult = XGBClf.predict(xFeats)
		# print type(_bdtResult)
		# print np.shape(_bdtResult)
		# print np.size(_bdtResult)
		bdtResult[0] = _bdtResult[0]

		_bdtProb = XGBClf.predict_proba(xFeats)
		sigProb[0] = _bdtProb.item(1)

		bgProb[0] = _bdtProb.item(0)

		# print 'sig: %s, bg: %s'%(sigProb[0], bgProb[0])

		_iEntry[0] = iEntry

		outTree.Fill()

	
	outTree.BuildIndex("iEntry")
	outFile.Write()
	outFile.Close()

	inputFile.Close()

	print 'BDT results saved in friend tree '+outFilePath


def main():
	now = datetime.datetime.now()
	print str(now)
	print 'This program will evaluate XGBClassifier predictions on TTrees and write results in friend trees.'
	print 'Model file : ' + model_file
	
	with open(sample_list, "r") as f:
		samples = [line.strip() for line in f if line.strip()]

	print 'Sample list: '
	for line in samples: print '\t' + line


	try:
		os.makedirs(output_directory)
	except OSError:
		if not os.path.isdir(output_directory):
			raise

	print 'Friend trees will be saved in ' + output_directory

	# load model
	clf = pk.load(open(model_file, 'rb'))
	# clf = XGBClassifier()
	# booster = Booster()
	# booster.load_model(model_file)
	# clf._Booster = booster
	# clf._le = LabelEncoder().fit([0, 1])

	for line in samples:
	  runBDTonTree(line, clf, output_directory)

	print 'Complete!'

main()
