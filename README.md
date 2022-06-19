# HLTTutorial

## Setup
Setup the release, import HLTrigger package, compile
```
cmsrel CMSSW_12_3_0 
cd CMSSW_12_3_0/src
cmsenv
git cms-addpkg HLTrigger/Configuration
git fetch official-cmssw pull/37552/head:pull37552
git checkout pull37552 HLTrigger/Configuration/python/customizeHLTforCMSSW.py
scram b -j 12
```
Creat a proxy (to access remote files) 

```
voms-proxy-init --voms cms
```

Clone this repository, and compile (N.B. some unused variables in the code introduce warnings)

```
git clone https://github.com/pallabidas/HLTTutorial.git -b 12_3_X
scram b -j 12
```
Move to HLTrigger/Configuration/test

```
cd HLTrigger/Configuration/test
```

Now, fetch the needed trigger configuration associated to the B-tagged Jet + soft muon path: 

```
hltGetConfiguration --cff --offline /dev/CMSSW_12_3_0/GRun --paths HLTriggerFirstPath,HLTriggerFinalPath --unprescale --l1 L1Menu_Collisions2018_v2_1_0-d1_xml --l1-emulator > HLT_TutoEffcySession_cff.py

hltGetConfiguration --cff /users/pdas/Jet110/HLT/V1 --globaltag auto:run3_hlt_GRun --unprescale --l1 L1Menu_Collisions2018_v2_1_0-d1_xml --l1-emulator >> HLT_TutoEffcySession_cff.py
 ```
 You then need to modify HLT_TutoEffcySession_cff.py : search for these two lines, they appear twice in the config file. Comment out their second appearance (not the first one). 
 ```
#/users/pdas/Jet110/HLT/V1 (CMSSW_12_3_0)   (already commented)                                                                                                                                       
#import FWCore.ParameterSet.Config as cms                                                                                                                                                                   
#fragment = cms.ProcessFragment( "HLT" )  
 ```
  
Copy the files to the python subdirectory
```
cp HLT_TutoEffcySession_cff.py ../python/.
```

 Produce a HLT2_HLT.py config file to be run with cmsRun HLT2_HLT.py. We will need to rerun HLT using RAW data and and then access offline information from MINIAOD. In order to do that, we need to specify the MINIAOD file as the main input file and the list of all its parent RAW files as secondary input files: 
```
cmsDriver.py HLT2 --step=HLT:TutoEffcySession --era=Run3 --mc --conditions auto:upgrade2021 --filein root://cms-xrd-global.cern.ch://store/mc/Run3Summer21MiniAOD/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v5-v1/30000/2ad0e924-f47b-4c16-89d3-b02164f10fc0.root --secondfilein root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/0a2412e5-e027-4fb0-b8fe-6faa36739bcc.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/cdf465f1-be67-4549-88f0-c3f526e447e3.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30005/c7dec4fe-8bfc-4332-bade-b0c43708af16.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/386da34a-36b8-4d1c-b1b5-f9c0e3e469ad.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/61e658ea-3786-4a36-bd10-6ab49b8bb3ce.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/7113b506-f65f-44ba-a8b4-1f36b5e6036d.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/ee1ff7f5-1bff-47b3-848d-84e0f73078d0.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30002/fd1c5d7b-fd14-4e12-83f3-a731090fc752.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/8fcd37dd-973e-4d83-bcda-a666178399da.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30004/f3916392-fc8c-44d8-ab44-940ad0e585cd.root,root://cms-xrd-global.cern.ch://store/mc/Run3Summer21DR/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW/FlatPU0to80FEVT_castor_120X_mcRun3_2021_realistic_v6-v1/30005/4e0dc09e-bfad-4179-a96e-b968b7e54c4c.root --processName=HLT2 -n 100 --no_exec
```

So far the HLT2_HLT.py file contains a configuration to rerun HLT on top of RAW. 
We need to add our EDAnalyzer as a second module in the process. After: 
```
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
```
Add: 
```
process.demo = cms.EDAnalyzer('TriggerAnalyzerRAWMiniAOD')
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( "out.root" )
                                   )
process.demo_step = cms.EndPath(process.demo)
```
Replace the line:
```
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])
```
by:
```
process.schedule.extend([process.endjob_step, process.demo_step])
```
N.B.: Keeping RECOSIMoutput_step in the process would create an updated (big!) RAW file 
with the information related to the trigger menu you reran. 
This file is typically quite big so we will drop it here. As an exercise, you may wish to produce 
it and take a look at its content, for example by doing: 
```
edmDumpEventContent --regex HLT HLT2_HLT.root
```
Test the code: 
```
cmsRun HLT2_HLT.py 
```

