#!/usr/bin/env python

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys

config = config()


#**************************submit function***********************
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
def submit(config):
	try:
		crabCommand('submit', config = config)
	except HTTPException as hte:
		print "Failed submitting task: %s" % (hte.headers)
	except ClientException as cle:
		print "Failed submitting task: %s" % (cle)
#****************************************************************


workarea='/afs/hep.wisc.edu/home/rpeterse/private/CMSSW_9_4_16_UL/src/ggAnalysis/ggNtuplizer/test/crab_submit//SinglePhotonRun2017E17Nov2017v1AOD/'

mainOutputDir = '/store/user/rpeterse/ggNtuplizer_METv2_AOD'


config.General.requestName = 'SinglePhotonRun2017E17Nov2017v1AOD'
config.General.transferLogs = True
config.General.workArea = '%s' % workarea


config.Site.storageSite = 'T2_US_Wisconsin'
config.Site.whitelist = ['T3_US_UCR','T3_US_FNALLPC','T2_US_Purdue','T3_US_Rice','T3_US_Rutgers','T3_US_FIT','T3_US_PSC','T3_US_OSU','T3_US_TAMU','T3_US_UMD','T3_US_VC3_NotreDame','T3_US_SDSC','T3_US_Colorado','T3_US_OSG','T3_US_Princeton_ICSE','T3_US_NERSC','T3_US_Baylor','T2_US_Nebraska','T2_US_UCSD','T2_US_Wisconsin','T2_US_MIT','T3_US_TACC','T3_US_TTU','T3_US_UMiss','T2_US_Caltech']
config.Site.blacklist = ['T2_US_Florida','T2_US_Vanderbilt','T3_US_PuertoRico']


config.JobType.psetName  = '/afs/hep.wisc.edu/home/rpeterse/private/CMSSW_9_4_16_UL/src/ggAnalysis/ggNtuplizer/test/run_data2017AOD_94X.py'
config.JobType.pluginName  = 'Analysis'


config.Data.inputDBS = 'global'
config.Data.inputDataset = '/SinglePhoton/Run2017E-17Nov2017-v1/AOD'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '%s' % mainOutputDir
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 2
config.Data.ignoreLocality = True
#config.Data.totalUnits = #totalUnits
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
submit(config)
