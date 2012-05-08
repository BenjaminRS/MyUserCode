import FWCore.ParameterSet.Config as cms

process = cms.Process("TrigMatch")

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
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	'/store/mc/Fall11/WH115_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/0227FFAF-7065-E111-AFFC-0030486361BC.root',
	'/store/mc/Fall11/WH115_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/50104CB2-3065-E111-B589-001A64789DEC.root',
	'/store/mc/Fall11/WH115_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/1E86BC3C-0665-E111-A2B8-001A647894F8.root',
	'/store/mc/Fall11/WH115_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/0C9DDCA1-3A65-E111-8874-001A647894F8.root'
#
#	'/store/mc/Fall11/ZH600_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/1C79B2EC-6465-E111-9E92-003048D45FDA.root'
#	'/store/mc/Fall11/ZH600_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/222352D9-D364-E111-B5E7-003048D45FC4.root'
#	'/store/mc/Fall11/ZH600_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/168CA83C-A365-E111-BFD3-0025B3E06468.root'
#	'/store/mc/Fall11/ZH600_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/1C9CD13E-4465-E111-8AB9-003048D47742.root'
    )
)

process.TrigMatch = cms.EDAnalyzer('TrigMatch')

#import HLTrigger.HLTfilters.hltHighLevel_cfi
#process.skimHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.skimHLTFilter.HLTPaths = cms.vstring("HLT_*")
#process.skimHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.TrigMatch.trigEventTag = cms.InputTag("hltTriggerSummaryAOD","","HLT")
#process.TrigMatch.filterName = cms.string("hltL1sL1DoubleEG125")
#process.TrigMatch.filterName = cms.string("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter"
#process.TrigMatch.filterName = cms.string("hltSingleMu40L2QualL3Filtered40")
process.TrigMatch.filterName = cms.string("hltSingleMu30L3Filtered30")
#process.TrigMatch.filterName = cms.string("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter")
#process.TrigMatch.filterName = cms.string("hltEG17EtFilter")
#process.TrigMatch.pathName = cms.string("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9")
process.TrigMatch.pathName = cms.string("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2")
#process.TrigMatch.pathName = cms.string("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")
#process.TrigMatch.pathName = cms.string("HLT_IsoMu30_v1")

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('LeptonTrigObjsEff-WH115-270412.root')
)

process.p = cms.Path(process.TrigMatch)
