###############################################
####
####   Jet Substructure Toolbox (jetToolBox)
####   Python function for easy access of jet substructure tools implemented in CMS
####
####   Alejandro Gomez Espinosa (gomez@physics.rutgers.edu)
####   First version: 2015 - 01 - 28
####
###############################################
import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from CommonTools.ParticleFlow.pfNoPileUpJME_cff  import *


def jetToolbox( proc, jetType, jetSequence, outputFile, minPt=100., 
		addPruning=False, zCut=0.1, rCut=0.5, addPrunedSubjets=False,
		addTrimming=False, rFiltTrim=0.2, ptFrac=0.03,
		addFiltering=False, rfilt=0.3, nfilt=3,
		addCMSTopTagger=False,
		addMassDrop=False,
		addHEPTopTagger=False,
		addNsub=False, addNsubUpTo5=False, 
		addQJets=False, 
		addSubjets=False, 
		miniAOD=True ):
	
	###############################################################################
	#######  Just defining simple variables
	###############################################################################
	supportedJetAlgos = { 'ak': 'AntiKt', 'ca' : 'CambridgeAachen', 'kt' : 'Kt' }
	recommendedJetAlgos = [ 'ak4', 'ak8', 'ca4', 'ca8' ]
	jetAlgo = ''
	algorithm = ''
	size = ''
	for type, tmpAlgo in supportedJetAlgos.iteritems(): 
		if type in jetType.lower():
			jetAlgo = type
			algorithm = tmpAlgo
			size = jetType.replace( type, '' )
	if jetAlgo == '': print 'Unsupported jet algorithm. Please use something like: jetType = CA8'

	jetSize = 0.
	if int(size) in range(0, 20): jetSize = int(size)/10.
	else: print 'jetSize has not a valid value. Insert a number between 1 and 20 after algorithm, like: AK8'
	### Trick for uppercase/lowercase algo name
	jetALGO = jetAlgo.upper()+size
	jetalgo = jetAlgo.lower()+size
	if( int(size) > 10 ): size = '10' 	### For JEC for jets larger than 1 
	recommended=False
	if jetalgo not in recommendedJetAlgos : print 'CMS recommends the following jet algoritms:', recommendedJetAlgos, 'You are using ', jetalgo
	else: recommended=True


	#################################################################################
	####### Toolbox start 
	#################################################################################

	elemToKeep = []
	jetSeq = cms.Sequence()
	genParticlesLabel = ''
	pvLabel = ''
	tvLabel = ''

	#### For MiniAOD
	if miniAOD:

		print '-------------- JETTOOLBOX RUNNING ON MiniAOD  ------------------'

		genParticlesLabel = 'prunedGenParticles'
		pvLabel = 'offlineSlimmedPrimaryVertices'
		tvLabel = 'unpackedTracksAndVertices'

		setattr( proc, 'chs', cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV')) )
		jetSeq += getattr(proc, 'chs')

		setattr( proc, jetalgo+'PFJetsCHS', 
				ak4PFJetsCHS.clone( src = 'chs', 
					doAreaFastjet = True, 
					rParam = jetSize, 
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )
		#elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHS_*_*' ]

		## Filter out neutrinos from packed GenParticles
		setattr( proc, 'packedGenParticlesForJetsNoNu', 
				cms.EDFilter("CandPtrSelector", 
					src = cms.InputTag("packedGenParticles"), 
					cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
					))
		jetSeq += getattr(proc, 'packedGenParticlesForJetsNoNu' )
		    
		setattr( proc, jetalgo+'GenJets', 
				ak4GenJets.clone( src = 'packedGenParticlesForJetsNoNu', 
					rParam = jetSize, 
					jetAlgorithm = algorithm ) ) 
		jetSeq += getattr(proc, jetalgo+'GenJets' )
		fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

	#### For AOD
	else:
		print '-------------- JETTOOLBOX RUNNING ON AOD  ------------------'

		proc.load('RecoJets.Configuration.GenJetParticles_cff')
		proc.load('RecoJets.Configuration.RecoPFJets_cff')
		setattr( proc, jetalgo+'GenJets', ak4GenJets.clone( src = 'genParticlesForJetsNoNu', rParam = jetSize, jetAlgorithm = algorithm ) ) 
		jetSeq += getattr(proc, jetalgo+'GenJets' )
		if not recommended: setattr( proc, jetalgo+'PFJetsCHS', ak4PFJets.clone( rParam = jetSize, jetAlgorithm = algorithm ) ) 
		jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )
		
		genParticlesLabel = 'genParticles'
		pvLabel = 'offlinePrimaryVertices'
		tvLabel = 'generalTracks'
		


	####  Creating PATjets
	proc.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
	#for Inclusive Vertex Finder
	proc.load("RecoBTag/Configuration/RecoBTag_cff")
	proc.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
	proc.inclusiveVertexFinder.tracks = cms.InputTag("unpackedTracksAndVertices")
	proc.inclusiveVertexFinder.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
	proc.trackVertexArbitrator.tracks = cms.InputTag("unpackedTracksAndVertices")
	proc.trackVertexArbitrator.primaryVertices = cms.InputTag("unpackedTracksAndVertices")

	#new input for impactParameterTagInfos, softleptons, IVF
	proc.impactParameterTagInfos.jetTracks = cms.InputTag("jetTracksAssociatorAtVertexSlimmedJetsAK8BTagged")
	proc.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
	proc.inclusiveVertexFinder.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
	proc.trackVertexArbitrator.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
	proc.softPFMuonsTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
	proc.softPFElectronsTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
	proc.softPFMuonsTagInfos.jets = cms.InputTag("patJetsSlimmedJetsAK8BTagged")
	proc.softPFElectronsTagInfos.jets = cms.InputTag("patJetsSlimmedJetsAK8BTagged") 
	proc.inclusiveSecondaryVertexFinderTagInfosV2 = proc.inclusiveSecondaryVertexFinderTagInfos.clone()
	proc.inclusiveSecondaryVertexFinderTagInfosV2.trackSelection.qualityClass = cms.string('any')

	## b-tag discriminators
	bTagDiscriminators = [
	    'trackCountingHighEffBJetTags',
	    'trackCountingHighPurBJetTags',
	    'jetProbabilityBJetTags',
	    'jetBProbabilityBJetTags',
	    'simpleSecondaryVertexHighEffBJetTags',
	    'simpleSecondaryVertexHighPurBJetTags',
	    'combinedSecondaryVertexBJetTags',
	    #'combinedInclusiveSecondaryVertexV2BJetTags'
	    ]

	addJetCollection(
			proc,
			labelName = jetALGO+'PFCHS',
			jetSource = cms.InputTag( jetalgo+'PFJetsCHS'),
			algo = jetalgo,
			rParam = jetSize,
			jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			trackSource = cms.InputTag( tvLabel ), 
			pvSource = cms.InputTag( pvLabel ), #'offlineSlimmedPrimaryVertices'),
			btagDiscriminators = bTagDiscriminators,
			getJetMCFlavour = False,
			outputModules = ['outputFile']
			) 

	#getattr( proc, 'patJets'+jetALGO+'PFCHS' ).addJetCharge = False 
	#getattr( proc, 'patJets'+jetALGO+'PFCHS' ).addAssociatedTracks = False 
	getattr( proc, 'patJetGenJetMatch'+jetALGO+'PFCHS' ).matched = cms.InputTag( jetalgo+'GenJets' ) 
	getattr( proc, 'patJetPartonMatch'+jetALGO+'PFCHS' ).matched = cms.InputTag( genParticlesLabel )  # 'prunedGenParticles' 
	getattr( proc, 'patJetCorrFactors'+jetALGO+'PFCHS' ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
	getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PFCHS' ).tracks = tvLabel  # 'unpackedTracksAndVertices'
	elemToKeep += [ 'keep *_patJets'+jetALGO+'PFCHS_*_*' ]
	jetSeq += getattr(proc, 'patJetGenJetMatch'+jetALGO+'PFCHS' )
	jetSeq += getattr(proc, 'patJetPartonMatch'+jetALGO+'PFCHS' )
	jetSeq += getattr(proc, 'patJetCorrFactors'+jetALGO+'PFCHS' )


	if addPruning: 

		if miniAOD: 
			setattr( proc, jetalgo+'PFJetsCHSPruned', 
				ak8PFJetsCHSPruned.clone( src = 'chs', 
					rParam = jetSize, 
					jetAlgorithm = algorithm, 
					zcut=zCut, 
					rcut_factor=rCut,
					jetCollInstanceName = 'Subjets') )
			setattr( proc, jetalgo+'PFJetsCHSPrunedLinks', 
				ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
					matched = cms.InputTag( jetalgo+'PFJetsCHSPruned'), 
					distMax = cms.double( jetSize ) ) )
		else:
			if not recommended:
				setattr( proc, jetalgo+'PFJetsCHSPruned', 
					ak8PFJetsCHSPruned.clone( 
						rParam = jetSize, 
						jetAlgorithm = algorithm, 
						zcut=zCut, 
						rcut_factor=rCut,
						jetCollInstanceName = 'Subjets') )
				setattr( proc, jetalgo+'PFJetsCHSPrunedLinks', 
					ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
						matched = cms.InputTag( jetalgo+'PFJetsCHSPruned'), 
						distMax = cms.double( jetSize ) ) )
			else: pass

		elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHSPrunedLinks_*_*'] 
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSPruned' )
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSPrunedLinks' )
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ jetalgo+'PFJetsCHSPrunedLinks']

		if addPrunedSubjets:
			setattr( proc, jetalgo+'GenJetsNoNuPruned',
					ak4GenJets.clone(
						SubJetParameters,
						rParam = jetSize,
						src = 'packedGenParticlesForJetsNoNu',
						usePruning = cms.bool(True),
						writeCompound = cms.bool(True),
						jetCollInstanceName=cms.string('SubJets')
						))

			addJetCollection(
					proc,
					labelName = jetALGO+'PFCHSPrunedSubjets',
					jetSource = cms.InputTag( jetalgo+'PFJetsCHSPruned', 'Subjets'),
					algo = jetalgo,
					rParam = jetSize,
					jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
					pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
					btagDiscriminators = bTagDiscriminators,
					genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNuPruned','SubJets'),
					getJetMCFlavour = False,
					#outputModules = ['outputFile']
					) 
			#elemToKeep += [ 'keep *_patJets'+jetALGO+'PFCHSPrunedSubjets_*_*' ]
			#if hasattr( proc, 'jetTracksAssociatorAtVertex' + jetALGO +'PFCHSPrunedSubjets' ): 
			#	proc.jetTracksAssociatorAtVertexAK8PFCHSPrunedSubjets.tracks = cms.InputTag("unpackedTracksAndVertices")

			getattr( proc,'patJetPartonMatch'+jetALGO+'PFCHSPrunedSubjets').matched = cms.InputTag('prunedGenParticles')
			getattr( proc,'patJetCorrFactors'+jetALGO+'PFCHSPrunedSubjets' ).primaryVertices = 'offlineSlimmedPrimaryVertices' 
			getattr( proc,'patJets'+jetALGO+'PFCHSPrunedSubjets').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
			getattr( proc,'patJets'+jetALGO+'PFCHSPrunedSubjets').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD

			## Establish references between PATified fat jets and subjets using the BoostedJetMerger
			setattr( proc, 'selectedPatJets'+jetALGO+'PFCHSPrunedPacked', 
					cms.EDProducer("BoostedJetMerger",
						jetSrc=cms.InputTag('patJets'+jetALGO+'PFCHS'),
						subjetSrc=cms.InputTag('selectedPatJets'+jetALGO+'PFCHSPrunedSubjets')
						))
			elemToKeep += [ 'keep *_selectedPatJets'+jetALGO+'PFCHSPrunedPacked_*_*' ]


	if addTrimming:

		if miniAOD: 
			setattr( proc, jetalgo+'PFJetsCHSTrimmed', 
					ak8PFJetsCHSTrimmed.clone( src = 'chs',
						rParam = jetSize, 
						jetAlgorithm = algorithm,
						rFilt= rFiltTrim,
						trimPtFracMin= ptFrac) ) 
			setattr( proc, jetalgo+'PFJetsCHSTrimmedLinks', 
					ak8PFJetsCHSTrimmedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
						matched = cms.InputTag( jetalgo+'PFJetsCHSTrimmed'), 
						distMax = cms.double( jetSize ) ) )
		else:
			if not recommended:
				setattr( proc, jetalgo+'PFJetsCHSTrimmed', 
						ak8PFJetsCHSTrimmed.clone( src = 'chs',
							rParam = jetSize, 
							jetAlgorithm = algorithm,
							rFilt= rFiltTrim,
							trimPtFracMin= ptFrac) ) 
				setattr( proc, jetalgo+'PFJetsCHSTrimmedLinks', 
						ak8PFJetsCHSTrimmedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
							matched = cms.InputTag( jetalgo+'PFJetsCHSTrimmed'), 
							distMax = cms.double( jetSize ) ) )
			else: pass

		elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHSTrimmedLinks_*_*'] 
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSTrimmed' )
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSTrimmedLinks' )
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ jetalgo+'PFJetsCHSTrimmedLinks']

	if addFiltering:

		setattr( proc, jetalgo+'PFJetsCHSFiltered', 
				ak8PFJetsCHSFiltered.clone( #src = 'chs', 
					rParam = jetSize, 
					jetAlgorithm = algorithm,
					rFilt= rfilt,
					nFilt= nfilt ) ) 
		if miniAOD: getattr( proc, jetalgo+'PFJetsCHSFiltered').src = 'chs'
		setattr( proc, jetalgo+'PFJetsCHSFilteredLinks', 
				ak8PFJetsCHSFilteredLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
					matched = cms.InputTag( jetalgo+'PFJetsCHSFiltered'), 
					distMax = cms.double( jetSize ) ) )
		elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHSFilteredLinks_*_*'] 
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSFiltered' )
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSFilteredLinks' )
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ jetalgo+'PFJetsCHSFilteredLinks']

	if addCMSTopTagger:

		if miniAOD: setattr( proc, 'cmsTopTagPFJetsCHS',  cmsTopTagPFJetsCHS.clone( src = 'chs' ) ) #, rParam = jetSize ) )
		setattr( proc, 'cmsTopTagPFJetsCHSLinks'+jetALGO, 
				ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
					matched = cms.InputTag("cmsTopTagPFJetsCHS"), 
					distMax = cms.double( jetSize ) ) )

		elemToKeep += [ 'keep *_cmsTopTagPFJetsCHSLinks'+jetALGO+'_*_*' ]
		jetSeq += getattr(proc, 'cmsTopTagPFJetsCHS' )
		jetSeq += getattr(proc, 'cmsTopTagPFJetsCHSLinks'+jetALGO )
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ 'cmsTopTagPFJetsCHSLinks'+jetALGO ]

	if ( 'CA' in jetALGO ) and addMassDrop :

		setattr( proc, jetalgo+'PFJetsCHSMassDropFiltered', ca15PFJetsCHSMassDropFiltered.clone( rParam = jetSize ) )
		if miniAOD: getattr( proc, jetalgo+'PFJetsCHSMassDropFiltered').src = 'chs'
		setattr( proc, jetalgo+'PFJetsCHSMassDropFilteredLinks', ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
			matched = cms.InputTag(jetalgo+'PFJetsCHSMassDropFiltered'), distMax = cms.double( jetSize ) ) )
		elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHSMassDropFilteredLinks_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ jetalgo+'PFJetsCHSMassDropFilteredLinks' ]
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSMassDropFiltered' )
		jetSeq += getattr(proc, jetalgo+'PFJetsCHSMassDropFilteredLinks' )
	else: 'CA is recommended for Mass Drop.'

	if ( 'CA' in jetALGO ) and ( jetSize > 1 ) and addHEPTopTagger: 

		if miniAOD: setattr( proc, 'hepTopTagPFJetsCHS', hepTopTagPFJetsCHS.clone( src = 'chs' ) )
		setattr( proc, 'hepTopTagPFJetsCHSLinks'+jetALGO, ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+"PFJetsCHS"), 
			matched = cms.InputTag("hepTopTagPFJetsCHS"), distMax = cms.double( jetSize ) ) )
		elemToKeep += [ 'keep *_hepTopTagPFJetsCHSLinks'+jetALGO+'_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ 'hepTopTagPFJetsCHSLinks'+jetALGO ]
		jetSeq += getattr(proc, 'hepTopTagPFJetsCHS' )
		jetSeq += getattr(proc, 'hepTopTagPFJetsCHSLinks'+jetALGO )
	else: 'CA and a jet Size bigger than 1 are recommended for HEPTopTagger.'

	####### Nsubjettiness
	if addNsub:
		from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

		if addNsubUpTo5: 
			setattr( proc, 'Njettiness'+jetALGO, 
					Njettiness.clone( src = cms.InputTag( jetalgo+'PFJetsCHS'), 
						cone = cms.double( jetSize ), 
						Njets = cms.vuint32(1,2,3,4,5) ) )
		else: 
			setattr( proc, 'Njettiness'+jetALGO, 
					Njettiness.clone( src = cms.InputTag( jetalgo+'PFJetsCHS'), 
						cone = cms.double( jetSize ) ) )

		elemToKeep += [ 'keep *_Njettiness'+jetALGO+'_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += ['Njettiness'+jetALGO+':tau1','Njettiness'+jetALGO+':tau2','Njettiness'+jetALGO+':tau3']  
		if addNsubUpTo5: getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += [ 'Njettiness'+jetALGO+':tau4','Njettiness'+jetALGO+':tau5' ]
		jetSeq += getattr(proc, 'Njettiness'+jetALGO )

	###### QJetsAdder
	if addQJets:
		### there must be a better way to do this random number introduction
		setattr( proc, 'RandomNumberGeneratorService', cms.Service("RandomNumberGeneratorService", 
							QJetsAdderCA8 = cms.PSet(initialSeed = cms.untracked.uint32(7)),
							QJetsAdderAK8 = cms.PSet(initialSeed = cms.untracked.uint32(31)),
							QJetsAdderCA15 = cms.PSet(initialSeed = cms.untracked.uint32(76)), ) )

		from RecoJets.JetProducers.qjetsadder_cfi import QJetsAdder
		setattr( proc, 'QJetsAdder'+jetALGO, 
				QJetsAdder.clone( src = cms.InputTag(jetalgo+'PFJetsCHS'), 
					jetRad = cms.double( jetSize ), 
					jetAlgo = cms.string( jetALGO[0:2] )))
		elemToKeep += [ 'keep *_QJetsAdder'+jetALGO+'_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PFCHS').userData.userFloats.src += ['QJetsAdder'+jetALGO+':QjetsVolatility']  
		jetSeq += getattr(proc, 'QJetsAdder'+jetALGO )

	####### Adding subjets
	jetSeq += getattr(proc, 'patJets'+jetALGO+'PFCHS' )
	if addSubjets : 
		setattr( proc, 'patJets'+jetALGO+'withSubjets', 
				cms.EDProducer('addSubjetProducer', 
					jets = cms.InputTag(jetalgo+'PFJetsCHS'), 
					patjets = cms.InputTag('patJets'+jetALGO+'PFCHS') ) )
		elemToKeep += [ 'keep *_patJets'+jetALGO+'withSubjets_*_*' ]
		jetSeq += getattr(proc, 'patJets'+jetALGO+'withSubjets' )
	

	### "return"
	setattr(proc, jetSequence, jetSeq)
	if hasattr(proc, outputFile): getattr(proc, outputFile).outputCommands += elemToKeep
	else: setattr( proc, outputFile, 
			cms.OutputModule('PoolOutputModule', 
				fileName = cms.untracked.string('jettoolbox.root'), 
				outputCommands = cms.untracked.vstring( elemToKeep ) ) )

