import FWCore.ParameterSet.Config as cms

process = cms.Process("BRS")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
#	fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/radburnsmith/Storage/044D42DB-A02B-E111-BB95-00261894388D.root'),
	fileNames = cms.untracked.vstring(),
        duplicateCheckMode=cms.untracked.string('checkEachFile')
)

process.brs = cms.EDAnalyzer('Ntuplize2',
		verbosity=cms.untracked.bool(False),
		TupFileName=cms.untracked.string("H2EJBRSNtuple-ZeeMC-v1.root"),
		RunningScheme=cms.untracked.string("ZEE_MC"),
		TriggerScheme=cms.untracked.string("electron")
)

process.p = cms.Path(process.brs)

