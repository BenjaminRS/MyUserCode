import FWCore.ParameterSet.Config as cms

process = cms.Process("BRS2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.INFO.limit = 0
process.MessageLogger.cout.threshold = cms.untracked.string('WARNING')
process.MessageLogger.cerr.FwkSummary = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(500),
    limit = cms.untracked.int32(10000000)
)       
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(500),
    limit = cms.untracked.int32(10000000)
)
import datetime
now = datetime.datetime.now()
import sys

numEvents=100
mass=115
step='A'

if len(sys.argv) < 5 or len(sys.argv) > 5:
	sys.exit('Usage: %s %s NumEvents Mass Step' % (sys.argv[0],sys.argv[1]))
elif len(sys.argv)==5:
	numEvents=int(sys.argv[2])
	mass=int(sys.argv[3])
	mass_str=sys.argv[3]
	step=sys.argv[4]
	print "Opening files with Mass =",mass,"and Step =",step,"then running over",numEvents,"events"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numEvents) )

from ZHFilesALL import *

process.source = cms.Source("PoolSource",
			    fileNames = GetEJFiles(mass,step),
#			fileNames = cms.untracked.vstring('dcap://heplnx209.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms/store/mc/Fall11/ZH115_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/0E99B607-A864-E111-8D51-00E08178C151.root'),
			    duplicateCheckMode=cms.untracked.string('checkEachFile')
)

process.brs = cms.EDAnalyzer('LJ')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('JetPtComp-EMF95Ak7_ZH'+mass_str+step+'_NoMetEJ_'+now.strftime("%d%m%y")+'.root'))

process.p = cms.Path(process.brs)
