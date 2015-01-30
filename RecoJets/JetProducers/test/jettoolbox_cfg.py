import FWCore.ParameterSet.Config as cms

process = cms.Process('jetToolbox')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

from RecoJets.JetProducers.jetToolbox_cff import *
#jetToolbox( process, 'ak6', 'ak6JetSubs', 'out', addSoftDrop=True ) #, addGroomers = False )
#jetToolbox( process, 'ca8', 'ca8JetSubs', 'out', addSubjets=True, addNsub=True, maxTau=6, addPruning=True, nfilt=5 , zCut=0.2, addTrimming=True, addCMSTopTagger=True, addHEPTopTagger=True, addMassDrop=True, addSoftDrop=True, addPUPPI=True, addCS=True, addSoftKiller=True )
#jetToolbox( process, 'ca15', 'ca15JetSubs', 'out' )
jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', addSubjets= True, miniAOD=False ) #, addGroomers = False )
jetToolbox( process, 'ca8', 'ca8JetSubs', 'out', miniAOD=False, addSubjets=True, addNsub=True, maxTau=6, addPruning=True, nfilt=5 , zCut=0.2 , addTrimming=True, addCMSTopTagger=True, addHEPTopTagger=True, addMassDrop=True, addSoftDrop=True, addPUPPI=True, addCS=True, addSoftKiller=True )
#jetToolbox( process, 'ca15', 'ca15JetSubs', 'out', miniAOD=False )


process.endpath = cms.EndPath(process.out)

process.source = cms.Source("PoolSource",
#		fileNames = cms.untracked.vstring('/store/user/jstupak/ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp/Spring14dr-PU_S14_POSTLS170_V6AN1-v1/140622_185946/0000/miniAOD-prod_PAT_1.root')
		fileNames = cms.untracked.vstring('/store/user/algomez/RPVSt100tojj_13TeV_pythia8_GENSIM/RPVSt100tojj_13TeV_pythia8_AODSIM_v706_PU40bx50/fffa6e37f9457f57c180027c6057d31f/RPVSt100tojj_13TeV_pythia8_AODSIM_PU40bx50_516_1_zGg.root')
#		fileNames = cms.untracked.vstring('/store/user/algomez/RPVSt100tojj_13TeV_pythia8_GENSIM/RPVSt100tojj_13TeV_pythia8_MiniAOD_v706_PU40bx25/b71e879835d2f0083a0e044b05216236/RPVSt100tojj_13TeV_pythia8_MiniAOD_PU40bx25_525_1_JVY.root')
		)
#
