import FWCore.ParameterSet.Config as cms
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

process = cms.Process("DarkShowerAnalysis")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("ScoutingParking.SctParAnalyzer.SctParLooper_cfi")
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun3_2024_realistic_v26')


#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '150X_mcRun3_2024_realistic_v2'


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-4p90-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2550000/9bb83424-f089-4474-906f-e0bc753c2863.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-4p90-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2550000/8e0880ea-9039-4c81-82ee-972c35e4670b.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-4p90-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2550000/36e3ed37-ce84-40c9-b583-f01f63b01152.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-4p90-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2550000/e5870e41-d269-48aa-af8d-bf10ab4bf939.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-1p90-mpi-4_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/140X_mcRun3_2024_realistic_v26-v2/2550000/65116876-069b-4907-8cf3-15631f41e7d3.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-1p90-mpi-4_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/140X_mcRun3_2024_realistic_v26-v2/2550000/da5ea65a-f251-442e-97f8-57c41e8143c5.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-1p90-mpi-4_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/249a7a55-8548-43af-a009-6dead32d9275.root',
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-1p90-mpi-4_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/140X_mcRun3_2024_realistic_v26-v2/2550000/c8b5db1d-59a0-4fc6-b7ec-253b84e7bf75.root'
        #'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-2p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/0cffd46a-92f5-4401-94d6-6ab55b6fd824.root',
        'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-2p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/4ef75fa6-c4d5-455b-871f-d2bd5a690dbc.root',
        'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-2p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/7bd6635a-8489-4e17-b161-fd83a7866aa2.root',
        'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-2p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/fddfed3a-699f-42b9-a94d-e6046af395ff.root',
        'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-2p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/c84603e5-2b4b-4420-be0f-569367249042.root',
        'root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAOD/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-2p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/2560000/917795b6-2b24-42f0-846f-3a351d38e829.root'
        
    )
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string("/afs/cern.ch/user/g/garciaja/private/output.root"))
#process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root"))

process.p = cms.Path(process.nTuplizer)
