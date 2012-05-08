## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# the source is already defined in patTemplate_cfg - overriding source and various other things
#process.load("CommonTools.ParticleFlow.Sources.source_ZtoEles_DBS_312_cfi")
process.source = cms.Source("PoolSource", 
        fileNames = cms.untracked.vstring(
#		'/store/data/Run2011A/SingleElectron/AOD/May10ReReco-v1/0005/F6AF587D-807B-E011-ADAD-0025901D4D6C.root'

#		'/store/mc/Fall11/ZH115_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/02EFCE31-3965-E111-B0A0-003048D47A62.root',
#		'/store/mc/Fall11/ZH125_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/0EEAAD6F-3365-E111-8B53-003048670B14.root',
#		'/store/mc/Fall11/ZH150_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/26BBBB31-5265-E111-A596-001A6478935C.root',
#		'/store/mc/Fall11/ZH200_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/18BDB898-4664-E111-A7E7-003048D46010.root',
#		'/store/mc/Fall11/ZH400_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/00F17FE1-8666-E111-BA39-003048D437CE.root',
#		'/store/mc/Fall11/ZH600_Dark_3step_25_12_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/08E69FB5-4665-E111-ABAB-00E08178C0C7.root',

#		'/store/mc/Fall11/ZH115_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/0695E903-5765-E111-8C7E-002590200AF4.root',
#		'/store/mc/Fall11/ZH125_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/108ACC41-8765-E111-9107-0025B3E05CF2.root',
#		'/store/mc/Fall11/ZH150_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/1C858925-5865-E111-AD2A-002590200978.root',
#		'/store/mc/Fall11/ZH200_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/02565D60-CA64-E111-927B-00E08178C17B.root',
#		'/store/mc/Fall11/ZH400_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/12ABCD3A-5165-E111-B65F-0025902009B4.root',
#		'/store/mc/Fall11/ZH600_Dark_3step_25_7_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/08F572A4-2A65-E111-A5E7-003048D46020.root',

#		'/store/mc/Fall11/ZH115_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/04A2E035-5B64-E111-AC1D-002481E14D64.root',
#		'/store/mc/Fall11/ZH125_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/1E874C65-9D64-E111-A97B-002481E14E2C.root',
#		'/store/mc/Fall11/ZH150_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/029A4EB4-6E64-E111-BD67-0025B3E05E1C.root',
#		'/store/mc/Fall11/ZH200_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/009B384F-9E65-E111-9CDB-003048D47A5E.root',
#		'/store/mc/Fall11/ZH400_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/106CE2F4-2C65-E111-A625-001A64789E1C.root',
#		'/store/mc/Fall11/ZH600_Dark_3step_15_6_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/06A20A23-8765-E111-A768-003048D45FE4.root',

		'/store/mc/Fall11/ZH115_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/04593EE6-9564-E111-8208-00E08179187F.root',
#		'/store/mc/Fall11/ZH125_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/027A5962-9464-E111-B6C6-003048635E2C.root',
#		'/store/mc/Fall11/ZH150_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/0C2EAEAB-6465-E111-B704-003048D46016.root',
#		'/store/mc/Fall11/ZH200_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/1AA0F1D2-5B65-E111-A71B-003048D45FF8.root',
#		'/store/mc/Fall11/ZH400_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/023B4203-4F64-E111-B58C-0030486361DC.root',
#		'/store/mc/Fall11/ZH600_Dark_3step_1_pt21_7TeV_Tune4C-pythia8/AODSIM/PU_S6_START44_V9B-v1/0000/02566683-4865-E111-81DA-003048D460F4.root',
	)
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.out.fileName = cms.untracked.string('SingleEle_7TeV_1kEvents_patTuple_PF2PAT_ak7-100412.root')
# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Configure PAT to use PF2PAT instead of AOD sources - this function will modify the PAT sequences. 
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK5"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=False, postfix=postfix) 

# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out.outputCommands = cms.untracked.vstring('drop *',
						   'keep *_genParticles*_*_*',
						   'keep *_*Jets*_*_*',
                                                   *patEventContentNoCleaning ) 

getattr(process,"patElectrons"+postfix).pfElectronSource = cms.InputTag("pfElectrons"+postfix)

# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False

# enable delta beta correction for muon selection in PF2PAT? 
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False

process.MessageLogger.cerr.FwkReport.reportEvery = 100


###################### Analyser:
process.brs = cms.EDAnalyzer("PatBasicAnalyzer",
	filename=cms.untracked.string("2604Jets-115D.csv"),
	mass=cms.untracked.uint32(115),
	step=cms.untracked.uint32(4)
#  electronSrc = cms.untracked.InputTag("selectedPatElectronsPFlow"),
#  jetSrc      = cms.untracked.InputTag("selectedPatJetsPFlow")
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Test_Jets2300412.root')
)
######################

process.p = cms.Path(
    getattr(process,"patPF2PATSequence"+postfix) + process.brs
)

