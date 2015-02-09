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
from PhysicsTools.PatAlgos.tools.jetTools import *
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from CommonTools.ParticleFlow.pfNoPileUpJME_cff  import *


def jetToolbox( proc, jetType, jetSequence, outputFile, PUMethod='CHS', miniAOD=True,
		minPt=100., 
		addPruning=False, zCut=0.1, rCut=0.5, addPrunedSubjets=False,
		addSoftDrop=False, betaCut=0.0,  zCutSD=0.1,
		addTrimming=False, rFiltTrim=0.2, ptFrac=0.03,
		addFiltering=False, rfilt=0.3, nfilt=3,
		addCMSTopTagger=False,
		addMassDrop=False,
		addHEPTopTagger=False,
		addNsub=False, maxTau=4, 
		addQJets=False, 
		addSubjets=False, 
		#addPUPPI=False, #sizePUPPI=0.4,
		#addCS=False,
		#addSoftKiller=False
		):
	
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
	if jetAlgo == '': print '|---- jetToolBox: Unsupported jet algorithm. Please use something like: jetType = CA8'

	jetSize = 0.
	if int(size) in range(0, 20): jetSize = int(size)/10.
	else: print '|---- jetToolBox: jetSize has not a valid value. Insert a number between 1 and 20 after algorithm, like: AK8'
	### Trick for uppercase/lowercase algo name
	jetALGO = jetAlgo.upper()+size
	jetalgo = jetAlgo.lower()+size

	### For JEC for jets larger than 1 
	if( int(size) > 10 ): 
		size = '10' 
		print '|--- For jets bigger than 1.0, the jet corrections are AK10PFchs.'

	recommended=False
	if jetalgo not in recommendedJetAlgos : print '|---- jetToolBox: CMS recommends the following jet algoritms:', recommendedJetAlgos, '. You are using', jetalgo,'.'
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

		print '|-------------- JETTOOLBOX RUNNING ON MiniAOD  ------------------'

		genParticlesLabel = 'prunedGenParticles'
		pvLabel = 'offlineSlimmedPrimaryVertices'
		svLabel = 'slimmedSecondaryVertices'
		tvLabel = 'unpackedTracksAndVertices'
		pfCand = 'packedPFCandidates'

		setattr( proc, 'chs', cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV')) )
		jetSeq += getattr(proc, 'chs')


		## Filter out neutrinos from packed GenParticles
		setattr( proc, 'packedGenParticlesForJetsNoNu', 
				cms.EDFilter("CandPtrSelector", 
					src = cms.InputTag("packedGenParticles"), 
					cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
					))
		jetSeq += getattr(proc, 'packedGenParticlesForJetsNoNu' )
		    
		setattr( proc, jetalgo+'GenJetsNoNu', 
				ak4GenJets.clone( src = 'packedGenParticlesForJetsNoNu', 
					rParam = jetSize, 
					jetAlgorithm = algorithm ) ) 
		jetSeq += getattr(proc, jetalgo+'GenJetsNoNu' )
		fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

		#for Inclusive Vertex Finder
		proc.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

	#### For AOD
	else:
		print '|-------------- JETTOOLBOX RUNNING ON AOD  ------------------'

		genParticlesLabel = 'genParticles'
		pvLabel = 'offlinePrimaryVertices'
		tvLabel = 'generalTracks'
		pfCand = 'particleFlow'
		svLabel = 'inclusiveSecondaryVertices'

		proc.load('RecoJets.Configuration.GenJetParticles_cff')
		proc.load('RecoJets.Configuration.RecoPFJets_cff')
		setattr( proc, jetalgo+'GenJetsNoNu', ak4GenJets.clone( src = 'genParticlesForJetsNoNu', rParam = jetSize, jetAlgorithm = algorithm ) ) 
		jetSeq += getattr(proc, jetalgo+'GenJetsNoNu' )
		if not recommended: setattr( proc, jetalgo+'PFJetsCHS', ak4PFJets.clone( rParam = jetSize, jetAlgorithm = algorithm ) ) 
		jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )
		
		

	## b-tag discriminators
	bTagDiscriminators = [
			'pfTrackCountingHighEffBJetTags',
			'pfTrackCountingHighPurBJetTags',
			'pfJetProbabilityBJetTags',
			'pfJetBProbabilityBJetTags',
			'pfSimpleSecondaryVertexHighEffBJetTags',
			'pfSimpleSecondaryVertexHighPurBJetTags',
			'pfCombinedSecondaryVertexBJetTags'
			#'pfCombinedInclusiveSecondaryVertexV2BJetTags'
	    ]

	####  Creating PATjets
	if 'Puppi' in PUMethod:

		proc.load('CommonTools.PileupAlgos.Puppi_cff')
		from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi
		setattr( proc, jetalgo+'PFJetsPuppi', 
				ak4PFJetsPuppi.clone( doAreaFastjet = True, 
					rParam = jetSize, 
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		if miniAOD:
			puppi.candName = cms.InputTag('packedPFCandidates')
			puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
		jetSeq += getattr(proc, 'puppi' )
		jetSeq += getattr(proc, jetalgo+'PFJetsPuppi' )

	elif 'CS' in PUMethod:

		from RecoJets.JetProducers.ak4PFJetsCS_cfi import ak4PFJetsCS
		#from RecoJets.JetProducers.ak4PFJetsCS_cfi import ak4PFJetsCS
		setattr( proc, jetalgo+'PFJetsCS', 
				ak4PFJetsCS.clone( doAreaFastjet = True, 
					csRParam = cms.double(jetSize),
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		if miniAOD: getattr( proc, jetalgo+'PFJetsCS').src = 'chs'
		jetSeq += getattr(proc, jetalgo+'PFJetsCS' )

		setattr( proc, jetalgo+'PFJetsCSConstituents', ak8PFJetsCSConstituents.clone( src = cms.InputTag(jetalgo+'PFJetsCS') ) )
		jetSeq += getattr(proc, jetalgo+'PFJetsCS' )

	elif 'SK' in PUMethod:

		proc.load('CommonTools.PileupAlgos.softKiller_cfi')
		from RecoJets.JetProducers.ak4PFJetsSK_cfi import ak4PFJetsSK
		setattr( proc, jetalgo+'PFJetsSK', 
				ak4PFJetsSK.clone( rParam = jetSize, 
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		if miniAOD: getattr( proc, 'softKiller' ).PFCandidates = cms.InputTag('packedPFCandidates')
		jetSeq += getattr(proc, 'softKiller' )
		jetSeq += getattr(proc, jetalgo+'PFJetsSK' )
	
	else: 
		setattr( proc, jetalgo+'PFJetsCHS', 
				ak4PFJetsCHS.clone( #src = 'chs', 
					doAreaFastjet = True, 
					rParam = jetSize, 
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		if miniAOD: getattr( proc, jetalgo+'PFJetsCHS').src = 'chs'
		jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )
		#elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHS_*_*' ]

	addJetCollection(
			proc,
			labelName = jetALGO+'PF'+PUMethod,
			jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod),
			algo = jetalgo,
			rParam = jetSize,
			jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			pfCandidates = cms.InputTag( pfCand ),  #'packedPFCandidates'),
			svSource = cms.InputTag( svLabel ),   #'slimmedSecondaryVertices'),
			genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNu'),
			pvSource = cms.InputTag( pvLabel ), #'offlineSlimmedPrimaryVertices'),
			btagDiscriminators = bTagDiscriminators,
			outputModules = ['outputFile']
			) 

	getattr( proc, 'patJetCorrFactors'+jetALGO+'PF'+PUMethod ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
	getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PF'+PUMethod ).tracks = tvLabel  # 'unpackedTracksAndVertices'

	getattr(proc,'patJetPartons').particles = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
	getattr(proc,'patJetPartonMatch'+jetALGO+'PF'+PUMethod).matched = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
	if hasattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod):
		    getattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod).extSVCollection = cms.InputTag( svLabel ) #'slimmedSecondaryVertices')
	if miniAOD:
		getattr(proc,'patJets'+jetALGO+'PF'+PUMethod).addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
		getattr(proc,'patJets'+jetALGO+'PF'+PUMethod).addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD
	elemToKeep += [ 'keep *_patJets'+jetALGO+'PF'+PUMethod+'_*_*' ]
	elemToKeep += [ 'drop *_patJets'+jetALGO+'PF'+PUMethod+'_calo*_*' ]
	elemToKeep += [ 'drop *_patJets'+jetALGO+'PF'+PUMethod+'_tagInfos_*' ]

	if addSoftDrop: 

		#if miniAOD or not recommended:
		setattr( proc, jetalgo+'PFJets'+PUMethod+'SoftDrop', 
			ak8PFJetsCHSSoftDrop.clone( 
				rParam = jetSize, 
				jetAlgorithm = algorithm, 
				useExplicitGhosts=True,
				zcut=zCutSD, 
				beta=betaCut,
				writeCompound = cms.bool(True),
				jetCollInstanceName=cms.string('SubJets') ) )
		if miniAOD: getattr( proc, jetalgo+'PFJets'+PUMethod+'SoftDrop').src = 'chs'
		setattr( proc, jetalgo+'PFJets'+PUMethod+'SoftDropLinks', 
			ak8PFJetsCHSSoftDropLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
				matched = cms.InputTag( jetalgo+'PFJets'+PUMethod+'SoftDrop'), 
				distMax = cms.double( jetSize ) ) )

		elemToKeep += [ 'keep *_'+jetalgo+'PFJets'+PUMethod+'SoftDropLinks_*_*'] 
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'SoftDrop' )
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'SoftDropLinks' )
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ jetalgo+'PFJets'+PUMethod+'SoftDropLinks']

	if addPruning or addPrunedSubjets: 

		#if miniAOD or not recommended:
		setattr( proc, jetalgo+'PFJets'+PUMethod+'Pruned', 
			ak8PFJetsCHSPruned.clone( 
				rParam = jetSize, 
				jetAlgorithm = algorithm, 
				zcut=zCut, 
				rcut_factor=rCut,
				writeCompound = cms.bool(True),
				doAreaFastjet = cms.bool(True),
				jetCollInstanceName=cms.string('SubJets') ) )
		if miniAOD: getattr( proc, jetalgo+'PFJets'+PUMethod+'Pruned').src = 'chs'
		setattr( proc, jetalgo+'PFJets'+PUMethod+'PrunedLinks', 
			ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
				matched = cms.InputTag( jetalgo+'PFJets'+PUMethod+'Pruned'), 
				distMax = cms.double( jetSize ) ) )

		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'Pruned' )
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'PrunedLinks' )
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ jetalgo+'PFJets'+PUMethod+'PrunedLinks']
		elemToKeep += [ 'keep *_'+jetalgo+'PFJets'+PUMethod+'PrunedLinks_*_*'] 

		if addPrunedSubjets:
			setattr( proc, jetalgo+'GenJetsNoNuPruned',
					ak4GenJets.clone(
						SubJetParameters,
						rParam = jetSize,
						usePruning = cms.bool(True),
						writeCompound = cms.bool(True),
						jetCollInstanceName=cms.string('SubJets')
						))
			if miniAOD: getattr( proc, jetalgo+'GenJetsNoNuPruned' ).src = 'packedGenParticlesForJetsNoNu'
			jetSeq += getattr(proc, jetalgo+'GenJetsNoNuPruned' )

			addJetCollection(
					proc,
					labelName = jetALGO+'PF'+PUMethod+'Pruned',
					jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod+'Pruned'),
					#algo = jetalgo,
					#rParam = jetSize,
					jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
					btagDiscriminators = ['None'],
					genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNu'),
					getJetMCFlavour = False,
					outputModules = ['outputFile']
					) 
			getattr( proc, 'patJetCorrFactors'+jetALGO+'PF'+PUMethod+'Pruned' ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
			getattr(proc,'patJetPartonMatch'+jetALGO+'PF'+PUMethod+'Pruned').matched = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
			#jetSeq += getattr(proc, jetALGO+'PF'+PUMethod+'Pruned' )

			addJetCollection(
					proc,
					labelName = jetALGO+'PF'+PUMethod+'PrunedSubjets',
					jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod+'Pruned', 'SubJets'),
					#algo = jetalgo,
					#rParam = jetSize,
					jetCorrections = ( 'AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
					pfCandidates = cms.InputTag( pfCand ),  #'packedPFCandidates'),
					pvSource = cms.InputTag( pvLabel), #'offlineSlimmedPrimaryVertices'),
					svSource = cms.InputTag( svLabel ),   #'slimmedSecondaryVertices'),
					btagDiscriminators = bTagDiscriminators,
					genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNuPruned','SubJets'),
					explicitJTA = True,
					getJetMCFlavour = False,
					svClustering = True,
					outputModules = ['outputFile']
					) 
			getattr( proc,'patJets'+jetALGO+'PF'+PUMethod+'PrunedSubjets').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
			getattr( proc,'patJets'+jetALGO+'PF'+PUMethod+'PrunedSubjets').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD
			if hasattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod+'PrunedSubjets'):
				    getattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod+'PrunedSubjets').extSVCollection = cms.InputTag( svLabel ) #'slimmedSecondaryVertices')
			getattr(proc,'patJetCorrFactors'+jetALGO+'PF'+PUMethod+'PrunedSubjets' ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
			getattr(proc,'patJetPartonMatch'+jetALGO+'PF'+PUMethod+'PrunedSubjets').matched = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')

			## Establish references between PATified fat jets and subjets using the BoostedJetMerger
			setattr( proc, 'patJets'+jetALGO+'PF'+PUMethod+'PrunedPacked', 
					cms.EDProducer("BoostedJetMerger",
						jetSrc=cms.InputTag('patJets'+jetALGO+'PF'+PUMethod+'Pruned'),
						subjetSrc=cms.InputTag('patJets'+jetALGO+'PF'+PUMethod+'PrunedSubjets')
						))
			jetSeq += getattr(proc, 'patJets'+jetALGO+'PF'+PUMethod+'PrunedPacked' )
			elemToKeep += [ 'keep *_patJets'+jetALGO+'PF'+PUMethod+'PrunedPacked_*_*' ]


	if addTrimming:

		#if miniAOD or not recommended:
		setattr( proc, jetalgo+'PFJets'+PUMethod+'Trimmed', 
				ak8PFJetsCHSTrimmed.clone( #src = 'chs',
					rParam = jetSize, 
					jetAlgorithm = algorithm,
					rFilt= rFiltTrim,
					trimPtFracMin= ptFrac) ) 
		if miniAOD: getattr( proc, jetalgo+'PFJets'+PUMethod+'Trimmed').src = 'chs'
		setattr( proc, jetalgo+'PFJets'+PUMethod+'TrimmedLinks', 
				ak8PFJetsCHSTrimmedLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
					matched = cms.InputTag( jetalgo+'PFJets'+PUMethod+'Trimmed'), 
					distMax = cms.double( jetSize ) ) )

		elemToKeep += [ 'keep *_'+jetalgo+'PFJets'+PUMethod+'TrimmedLinks_*_*'] 
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'Trimmed' )
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'TrimmedLinks' )
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ jetalgo+'PFJets'+PUMethod+'TrimmedLinks']

	if addFiltering:

		setattr( proc, jetalgo+'PFJets'+PUMethod+'Filtered', 
				ak8PFJetsCHSFiltered.clone( #src = 'chs', 
					rParam = jetSize, 
					jetAlgorithm = algorithm,
					rFilt= rfilt,
					nFilt= nfilt ) ) 
		if miniAOD: getattr( proc, jetalgo+'PFJets'+PUMethod+'Filtered').src = 'chs'
		setattr( proc, jetalgo+'PFJets'+PUMethod+'FilteredLinks', 
				ak8PFJetsCHSFilteredLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
					matched = cms.InputTag( jetalgo+'PFJets'+PUMethod+'Filtered'), 
					distMax = cms.double( jetSize ) ) )
		elemToKeep += [ 'keep *_'+jetalgo+'PFJets'+PUMethod+'FilteredLinks_*_*'] 
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'Filtered' )
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'FilteredLinks' )
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ jetalgo+'PFJets'+PUMethod+'FilteredLinks']

	if addCMSTopTagger:

		if miniAOD: setattr( proc, 'cmsTopTagPFJets'+PUMethod,  cmsTopTagPFJetsCHS.clone( src = 'chs' ) ) #, rParam = jetSize ) )
		setattr( proc, 'cmsTopTagPFJets'+PUMethod+'Links'+jetALGO, 
				ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
					matched = cms.InputTag('cmsTopTagPFJets'+PUMethod), 
					distMax = cms.double( jetSize ) ) )

		elemToKeep += [ 'keep *_cmsTopTagPFJets'+PUMethod+'Links'+jetALGO+'_*_*' ]
		jetSeq += getattr(proc, 'cmsTopTagPFJets'+PUMethod )
		jetSeq += getattr(proc, 'cmsTopTagPFJets'+PUMethod+'Links'+jetALGO )
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ 'cmsTopTagPFJets'+PUMethod+'Links'+jetALGO ]

	if ( 'CA' in jetALGO ) and addMassDrop :

		setattr( proc, jetalgo+'PFJets'+PUMethod+'MassDropFiltered', ca15PFJetsCHSMassDropFiltered.clone( rParam = jetSize ) )
		if miniAOD: getattr( proc, jetalgo+'PFJets'+PUMethod+'MassDropFiltered').src = 'chs'
		setattr( proc, jetalgo+'PFJets'+PUMethod+'MassDropFilteredLinks', ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
			matched = cms.InputTag(jetalgo+'PFJets'+PUMethod+'MassDropFiltered'), distMax = cms.double( jetSize ) ) )
		elemToKeep += [ 'keep *_'+jetalgo+'PFJets'+PUMethod+'MassDropFilteredLinks_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ jetalgo+'PFJets'+PUMethod+'MassDropFilteredLinks' ]
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'MassDropFiltered' )
		jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'MassDropFilteredLinks' )
	else: 'CA is recommended for Mass Drop.'

	if ( 'CA' in jetALGO ) and ( jetSize > 1 ) and addHEPTopTagger: 

		if miniAOD: setattr( proc, 'hepTopTagPFJets'+PUMethod, hepTopTagPFJetsCHS.clone( src = 'chs' ) )
		setattr( proc, 'hepTopTagPFJets'+PUMethod+'Links'+jetALGO, ak8PFJetsCHSPrunedLinks.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
			matched = cms.InputTag("hepTopTagPFJets'+PUMethod+'"), distMax = cms.double( jetSize ) ) )
		elemToKeep += [ 'keep *_hepTopTagPFJets'+PUMethod+'Links'+jetALGO+'_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [ 'hepTopTagPFJets'+PUMethod+'Links'+jetALGO ]
		jetSeq += getattr(proc, 'hepTopTagPFJets'+PUMethod )
		jetSeq += getattr(proc, 'hepTopTagPFJets'+PUMethod+'Links'+jetALGO )
	else: 'CA and a jet Size bigger than 1 are recommended for HEPTopTagger.'

	####### Nsubjettiness
	if addNsub:
		from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

		rangeTau = range(1,maxTau+1)
		setattr( proc, 'Njettiness'+jetALGO+PUMethod, 
				Njettiness.clone( src = cms.InputTag( jetalgo+'PFJets'+PUMethod), 
					Njets=cms.vuint32(rangeTau),         # compute 1-, 2-, 3-, 4- subjettiness
					# variables for measure definition : 
					measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
					beta = cms.double(1.0),              # CMS default is 1
					R0 = cms.double( jetSize ),              # CMS default is jet cone size
					Rcutoff = cms.double( -999.0),       # not used by default
					# variables for axes definition :
					axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
					nPass = cms.int32(-999),             # not used by default
					akAxesR0 = cms.double(-999.0) ) )        # not used by default

		elemToKeep += [ 'keep *_Njettiness'+jetALGO+PUMethod+'_*_*' ]
		for tau in rangeTau: getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += ['Njettiness'+jetALGO+PUMethod+':tau'+str(tau) ] 
		jetSeq += getattr(proc, 'Njettiness'+jetALGO+PUMethod )

	###### QJetsAdder
	if addQJets:
		#This is the old way before 731
		### there must be a better way to do this random number introduction
		setattr( proc, 'RandomNumberGeneratorService', cms.Service("RandomNumberGeneratorService", 
							QJetsAdderCA8 = cms.PSet(initialSeed = cms.untracked.uint32(7)),
							QJetsAdderAK8 = cms.PSet(initialSeed = cms.untracked.uint32(31)),
							QJetsAdderCA15 = cms.PSet(initialSeed = cms.untracked.uint32(76)), ) )

		from RecoJets.JetProducers.qjetsadder_cfi import QJetsAdder
		setattr( proc, 'QJetsAdder'+jetALGO, 
				QJetsAdder.clone( src = cms.InputTag(jetalgo+'PFJets'+PUMethod), 
					jetRad = cms.double( jetSize ), 
					jetAlgo = cms.string( jetALGO[0:2] )))
		elemToKeep += [ 'keep *_QJetsAdder'+jetALGO+'_*_*' ]
		getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += ['QJetsAdder'+jetALGO+':QjetsVolatility']  
		jetSeq += getattr(proc, 'QJetsAdder'+jetALGO )
		'''
		### This is for 731 or higher
		if 'ak4' in jetalgo:
			proc.load('RecoJets.JetProducers.QGTagger_cfi')
			proc.QGTagger.srcJets = cms.InputTag(jetalgo+'PFJets'+PUMethod)    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
			proc.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')        # Other options (might need to add an ESSource for it): see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
			elemToKeep += [ 'keep *_QGTagger_*_*' ]
			getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += ['QJetsAdder'+jetALGO+':qgLikelihood']  
			jetSeq += getattr(proc, 'QGTagger' )
		else:
			'QGTagger is optimized for ak4 jets.'
		'''

	'''
	####### Adding PUPI
	if addPUPPI:

		proc.load('CommonTools.PileupAlgos.Puppi_cff')
		from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi
		setattr( proc, jetalgo+'PFJetsPuppi', 
				ak4PFJetsPuppi.clone( doAreaFastjet = True, 
					rParam = jetSize, 
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		if miniAOD:
			puppi.candName = cms.InputTag('packedPFCandidates')
			puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
		jetSeq += getattr(proc, 'puppi' )
		jetSeq += getattr(proc, jetalgo+'PFJetsPuppi' )
		#elemToKeep += [ 'keep *_'+jetalgo+'PFJetsPuppi_*_*' ]

		addJetCollection(
				proc,
				labelName = jetALGO+'PFPuppi',
				jetSource = cms.InputTag( jetalgo+'PFJetsPuppi'),
				algo = jetalgo,
				rParam = jetSize,
				jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
				trackSource = cms.InputTag( tvLabel ), 
				pvSource = cms.InputTag( pvLabel ), #'offlineSlimmedPrimaryVertices'),
				btagDiscriminators = bTagDiscriminators,
				getJetMCFlavour = False,
				outputModules = ['outputFile']
				) 

		#getattr( proc, 'patJets'+jetALGO+'PFPuppi' ).addJetCharge = False 
		#getattr( proc, 'patJets'+jetALGO+'PFPuppi' ).addAssociatedTracks = False 
		getattr( proc, 'patJetGenJetMatch'+jetALGO+'PFPuppi' ).matched = cms.InputTag( jetalgo+'GenJetsNoNu' ) 
		getattr( proc, 'patJetPartonMatch'+jetALGO+'PFPuppi' ).matched = cms.InputTag( genParticlesLabel )  # 'prunedGenParticles' 
		getattr( proc, 'patJetCorrFactors'+jetALGO+'PFPuppi' ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
		getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PFPuppi' ).tracks = tvLabel  # 'unpackedTracksAndVertices'
		elemToKeep += [ 'keep *_patJets'+jetALGO+'PFPuppi_*_*' ]
		jetSeq += getattr(proc, 'patJetGenJetMatch'+jetALGO+'PFPuppi' )
		jetSeq += getattr(proc, 'patJetPartonMatch'+jetALGO+'PFPuppi' )
		jetSeq += getattr(proc, 'patJetCorrFactors'+jetALGO+'PFPuppi' )

	#### Adding Constituent Substraction
	if addCS:

		from RecoJets.JetProducers.ak4PFJetsCS_cfi import ak4PFJetsCS
		setattr( proc, jetalgo+'PFJetsCS', 
				ak4PFJetsCS.clone( doAreaFastjet = True, 
					csRParam = cms.double(jetSize),
					jetAlgorithm = algorithm,  
					jetPtMin = minPt )) 
		if miniAOD: getattr( proc, jetalgo+'PFJetsCS').src = 'chs'
		jetSeq += getattr(proc, jetalgo+'PFJetsCS' )

		addJetCollection(
				proc,
				labelName = jetALGO+'PFCS',
				jetSource = cms.InputTag( jetalgo+'PFJetsCS'),
				algo = jetalgo,
				rParam = jetSize,
				jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
				trackSource = cms.InputTag( tvLabel ), 
				pvSource = cms.InputTag( pvLabel ), #'offlineSlimmedPrimaryVertices'),
				btagDiscriminators = bTagDiscriminators,
				getJetMCFlavour = False,
				outputModules = ['outputFile']
				) 

		#getattr( proc, 'patJets'+jetALGO+'PFCS' ).addJetCharge = False 
		#getattr( proc, 'patJets'+jetALGO+'PFCS' ).addAssociatedTracks = False 
		getattr( proc, 'patJetGenJetMatch'+jetALGO+'PFCS' ).matched = cms.InputTag( jetalgo+'GenJetsNoNu' ) 
		getattr( proc, 'patJetPartonMatch'+jetALGO+'PFCS' ).matched = cms.InputTag( genParticlesLabel )  # 'prunedGenParticles' 
		getattr( proc, 'patJetCorrFactors'+jetALGO+'PFCS' ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
		getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PFCS' ).tracks = tvLabel  # 'unpackedTracksAndVertices'
		elemToKeep += [ 'keep *_patJets'+jetALGO+'PFCS_*_*' ]
		jetSeq += getattr(proc, 'patJetGenJetMatch'+jetALGO+'PFCS' )
		jetSeq += getattr(proc, 'patJetPartonMatch'+jetALGO+'PFCS' )
		jetSeq += getattr(proc, 'patJetCorrFactors'+jetALGO+'PFCS' )

	#### addind Soft Killer
	if addSoftKiller:

		if not miniAOD:
			proc.load('CommonTools.PileupAlgos.softKiller_cfi')
			from RecoJets.JetProducers.ak4PFJetsSK_cfi import ak4PFJetsSK
			setattr( proc, jetalgo+'PFJetsSK', 
					ak4PFJetsSK.clone( rParam = jetSize, 
						jetAlgorithm = algorithm,  
						jetPtMin = minPt )) 
			#if miniAOD: getattr( proc, 'softKiller' ).PFCandidates = cms.InputTag('packedPFCandidates')
			jetSeq += getattr(proc, 'softKiller' )
			jetSeq += getattr(proc, jetalgo+'PFJetsSK' )

			addJetCollection(
					proc,
					labelName = jetALGO+'PFSK',
					jetSource = cms.InputTag( jetalgo+'PFJetsSK'),
					algo = jetalgo,
					rParam = jetSize,
					jetCorrections = ( 'AK'+size+'PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
					trackSource = cms.InputTag( tvLabel ), 
					pvSource = cms.InputTag( pvLabel ), #'offlineSlimmedPrimaryVertices'),
					btagDiscriminators = bTagDiscriminators,
					getJetMCFlavour = False,
					outputModules = ['outputFile']
					) 

			#getattr( proc, 'patJets'+jetALGO+'PFSK' ).addJetCharge = False 
			#getattr( proc, 'patJets'+jetALGO+'PFSK' ).addAssociatedTracks = False 
			getattr( proc, 'patJetGenJetMatch'+jetALGO+'PFSK' ).matched = cms.InputTag( jetalgo+'GenJetsNoNu' ) 
			getattr( proc, 'patJetPartonMatch'+jetALGO+'PFSK' ).matched = cms.InputTag( genParticlesLabel )  # 'prunedGenParticles' 
			getattr( proc, 'patJetCorrFactors'+jetALGO+'PFSK' ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
			getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PFSK' ).tracks = tvLabel  # 'unpackedTracksAndVertices'
			elemToKeep += [ 'keep *_patJets'+jetALGO+'PFSK_*_*' ]
			jetSeq += getattr(proc, 'patJetGenJetMatch'+jetALGO+'PFSK' )
			jetSeq += getattr(proc, 'patJetPartonMatch'+jetALGO+'PFSK' )
			jetSeq += getattr(proc, 'patJetCorrFactors'+jetALGO+'PFSK' )

		else: print '|---- jetToolBox: SoftKiller needs PFCandidate collection, not available in miniAOD.'
	'''

	####### Adding subjets
	jetSeq += getattr(proc, 'patJets'+jetALGO+'PF'+PUMethod )
	if addSubjets : 
		setattr( proc, 'patJets'+jetALGO+PUMethod+'withSubjets', 
				cms.EDProducer('addSubjetProducer', 
					jets = cms.InputTag(jetalgo+'PFJets'+PUMethod), 
					patjets = cms.InputTag('patJets'+jetALGO+'PF'+PUMethod) ) )
		elemToKeep += [ 'keep *_patJets'+jetALGO+PUMethod+'withSubjets_*_*' ]
		jetSeq += getattr(proc, 'patJets'+jetALGO+PUMethod+'withSubjets' )
	

	### "return"
	setattr(proc, jetSequence, jetSeq)
	if hasattr(proc, outputFile): getattr(proc, outputFile).outputCommands += elemToKeep
	else: setattr( proc, outputFile, 
			cms.OutputModule('PoolOutputModule', 
				fileName = cms.untracked.string('jettoolbox.root'), 
				outputCommands = cms.untracked.vstring( elemToKeep ) ) )


