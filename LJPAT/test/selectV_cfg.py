import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
#    'file:/opt/ppd/scratch/radburnsmith/Storage/ZZFall11MC_patTuple_PF2PAT_TEST.root'
#	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/ZZFall11MC_patTuple_PF2PAT_TEST250112.root'
#	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/Reco110_patTuple_PF2PAT_TEST210212.root'
#	'file:/opt/ppd/scratch/radburnsmith/Storage/3stepC_15.0_6.0_reco_175_darkstuff.root'
#	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/3SAReco54_patTuple_PF2PAT_TEST020312.root'
#	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/3SBReco166_patTuple_PF2PAT_TEST020312.root'
#	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/3SCReco175_patTuple_PF2PAT_TEST290212.root'
#	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/3SDReco82_patTuple_PF2PAT_TEST020312.root'
	'file:/opt/ppd/newscratch/radburnsmith/CMSSW_4_4_2_patch9/src/LJets/LJPAT/ZH200SA_V9B_76A8_patTuple_PF2PAT_TEST080312.root'
  )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.brs = cms.EDAnalyzer("VSelector",
#  photonSrc   = cms.untracked.InputTag(""),
  electronSrc = cms.untracked.InputTag("selectedPatElectronsPFlow"),
#  electronSrc = cms.untracked.InputTag("selectedPatElectronsNonIsoPFlow"),
#  muonSrc     = cms.untracked.InputTag("selectedPatMuonsPFlow"),                                             
#  tauSrc      = cms.untracked.InputTag("selectedPatTausPFlow"),
#  jetSrc      = cms.untracked.InputTag("selectedPatJetsPFlow")
#  jetSrc      = cms.untracked.InputTag("selectedPatJetsPFlow2")
#  metSrc      = cms.untracked.InputTag("patMETsPFlow")
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('ElePts_130312.root')
)

process.p = cms.Path(process.brs)

