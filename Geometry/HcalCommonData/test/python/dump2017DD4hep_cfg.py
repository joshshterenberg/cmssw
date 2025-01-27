import FWCore.ParameterSet.Config as cms

process = cms.Process("DUMP")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 5
#if 'MessageLogger' in process.__dict__:
#   process.MessageLogger.PixelGeom=dict()
#   process.MessageLogger.TIBGeom=dict()
#   process.MessageLogger.TIDGeom=dict()
#   process.MessageLogger.TOBGeom=dict()
#   process.MessageLogger.TECGeom=dict()
#   process.MessageLogger.TGeoMgrFromDdd=dict()

process.DDDetectorESProducer = cms.ESSource("DDDetectorESProducer",
                                            confGeomXMLFiles = cms.FileInPath('Geometry/HcalCommonData/data/cmsExtendedGeometry2017.xml'),
                                            appendToDataLabel = cms.string('cms2017')
                                            )

process.testDump = cms.EDAnalyzer("DDTestDumpFile",
                                  outputFileName = cms.untracked.string('cms2017DD4hep.root'),
                                  DDDetector = cms.ESInputTag('','cms2017')
                                  )

process.p = cms.Path(process.testDump)
