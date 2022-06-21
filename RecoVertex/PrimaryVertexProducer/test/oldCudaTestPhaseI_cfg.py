import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Era_Run3_cff import Run3
#from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVerticesDumbFitter as offlinePrimaryVertices
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVerticesDumbFitter as offlinePrimaryVertices
#process = cms.Process("Demo", Run2_2018)
process = cms.Process("Demo", Run3)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')
process.load("HeterogeneousCore.CUDAServices.CUDAService_cfi")
process.load("HeterogeneousCore.CUDAServices.NVProfilerService_cfi")

#process.load( "HLTrigger.Timer.FastTimerService_cfi" )

#process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2023_realistic', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
"""
process.IgProfService = cms.Service("IgProfService",
  reportFirstEvent            = cms.untracked.int32(0),
  reportEventInterval         = cms.untracked.int32(25),
  reportToFileAtPostEvent     = cms.untracked.string("| gzip -c > igdqm.%I.gz")
)
"""
process.source = cms.Source("PoolSource",
# replace 'myfile.root' with the source file you want to use 
fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/user/g/gpizzati/CMSSW_12_2_0_pre3/src/step2.root' 
#'/store/relval/CMSSW_12_2_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_122X_mcRun4_realistic_v2_2026D88PU200-v1/2580000/aca7b050-5990-4576-a9ee-f41ac82e5b86.root'
#'/store/relval/CMSSW_12_1_0_pre5/RelValHSCPgluino_M-1000_TuneCP5_13TeV-pythia8/GEN-SIM-RECO/PU_121X_mcRun3_2021_realistic_v15_HighStat-v1/2580000/5073d52d-ece2-4bf3-ae32-baeec737b51d.root',
#"/store/relval/CMSSW_12_4_0_pre3/RelValTTbarToDilepton_14TeV/GEN-SIM-RECO/PU_123X_mcRun3_2021_realistic_v14-v1/2580000/482fcd10-ed9c-4163-b6b9-c5d45f9b38d7.root",
#"/store/relval/CMSSW_11_0_0_pre13/RelValTTbar_14TeV/GEN-SIM-RECO/PU_110X_mcRun3_2023_realistic_v6-v1/20000/BF639D16-A5A8-FB4C-8A8E-7E4A7CA41375.root"
#"/store/relval/CMSSW_12_0_0_pre4/RelValTTbar_14TeV/GEN-SIM/120X_mcRun3_2023_realistic_v2-v1/00000/65c9e814-dd71-48b4-9583-93009fc2f29d.root"
"/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun3_2021_realistic_v14-v1/2580000/16ca3cdd-33b1-457e-aacf-73732649acca.root"

#'/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/7781d089-b51a-495a-b1ba-384c15e90749.root'
#,'/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/f6b68ca4-5b0e-42bb-b1d0-f94480067693.root',
#'/store/relval/CMSSW_12_4_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v11_2026D88PU200-v1/2580000/876a46e3-477e-4c53-8a4a-c16e7c8dee0b.root'

#'/store/relval/CMSSW_12_3_0_pre5/RelValTTbar_14TeV/GEN-SIM-RECO/123X_mcRun4_realistic_v4_2026D88noPU-v1/10000/7cbdf153-c397-410e-8af5-dc3905565ced.root'
),
#firstEvent = cms.untracked.uint32(2)
skipEvents=cms.untracked.uint32(0)
#skipEvents=cms.untracked.uint32(87)
)
process.NVProfilerService = cms.Service("NVProfilerService"
)

"""
process.ThroughputService = cms.Service('ThroughputService',
    eventRange = cms.untracked.uint32(10),
    eventResolution = cms.untracked.uint32(1),
    printEventSummary = cms.untracked.bool(True),
    enableDQM = cms.untracked.bool(True),
    dqmPathByProcesses = cms.untracked.bool(False),
    dqmPath = cms.untracked.string('Throughput'),
    timeRange = cms.untracked.double(1000),
    timeResolution = cms.untracked.double(1)
)

process.MessageLogger.cerr.ThroughputService = cms.untracked.PSet(
    limit = cms.untracked.int32(10000000)
)
"""
"""
process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)
process.FastTimerService.jsonFileName = cms.untracked.string('resources.json')
"""
process.demo = offlinePrimaryVertices
"""

process.demo = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        algorithm   = cms.string("DA_vectCUDA"),
        TkDAClusParameters = cms.PSet(
            coolingFactor = cms.double(0.6),  # moderate annealing speed
            zrange = cms.double(4.),          # consider only clusters within 4 sigma*sqrt(T) of a track
            delta_highT = cms.double(1.e-2),  # convergence requirement at high T
            delta_lowT = cms.double(1.e-3),   # convergence requirement at low T
            convergence_mode = cms.int32(0),  # 0 = two steps, 1 = dynamic with sqrt(T)
            Tmin = cms.double(2.0),           # end of vertex splitting
            Tpurge = cms.double(2.0),         # cleaning 
            Tstop = cms.double(0.5),          # end of annealing 
            vertexSize = cms.double(0.006),   # added in quadrature to track-z resolutions        
            d0CutOff = cms.double(3.),        # downweight high IP tracks 
            dzCutOff = cms.double(3.),        # outlier rejection after freeze-out (T<Tmin)       
            zmerge = cms.double(1e-2),        # merge intermediat clusters separated by less than zmerge
            uniquetrkweight = cms.double(0.8),# require at least two tracks with this weight at T=Tpurge
            uniquetrkminp = cms.double(0.0)   # minimal a priori track weight for counting unique tracks
            )
    ),
    TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        maxNormalizedChi2 = cms.double(10.0),
        minPixelLayersWithHits=cms.int32(2),
        minSiliconLayersWithHits = cms.int32(5),
        maxD0Significance = cms.double(4.0), 
        maxD0Error = cms.double(1.0), 
        maxDzError = cms.double(1.0), 
        minPt = cms.double(0.0),
        maxEta = cms.double(2.4),
        trackQuality = cms.string("any")
    ),
    TrackLabel = cms.InputTag("generalTracks"), #TrackLabel = cms.InputTag("pixelTracks","","RECO"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    verbose = cms.untracked.bool(False),
    vertexCollections = cms.VPSet(
     [cms.PSet(label=cms.string(""),
               algorithm=cms.string("AdaptiveVertexFitter"),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(0.0),
               useBeamConstraint = cms.bool(False),
               maxDistanceToBeam = cms.double(1.0)
               ),
      cms.PSet(label=cms.string("WithBS"),
               algorithm = cms.string('AdaptiveVertexFitter'),
               chi2cutoff = cms.double(2.5),
               minNdof=cms.double(2.0),
               useBeamConstraint = cms.bool(True),
               maxDistanceToBeam = cms.double(1.0),
               )
      ]
    ),
#    onGPU = cms.bool(True),
    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag("")
)
"""
"""
process.demo = cms.EDProducer(
    'gpuTest',
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),

)
"""
process.out = cms.OutputModule("PoolOutputModule", 
    fileName= cms.untracked.string("file:test_gpu.root"),
    outputCommands = cms.untracked.vstring(
                                'drop *_*_*_*',
                                'keep *_demo_*_*'
    )
)

process.p = cms.Path(process.demo)
process.ep = cms.EndPath(process.out)
