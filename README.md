# HLTTutorial

## Setup
Setup the release, import HLTrigger package, compile
```
cmsrel CMSSW_11_2_0
cd CMSSW_11_2_0/src
cmsenv
git cms-merge-topic pallabidas:l1t-integration-test-17May
git cms-addpkg HLTrigger/Configuration
scram b -j 12
```
Creat a proxy (to access remote files) 

```
voms-proxy-init --voms cms --valid 168:00 
```

Clone this repository, and compile (N.B. some unused variables in the code introduce warnings)

```
git clone -b 2021_July1 git@github.com:pallabidas/HLTTutorial.git
scram b -j 12
```
Move to HLTrigger/Configuration/test

```
cd HLTrigger/Configuration/test
```

Now, fetch the needed trigger configuration associated to the B-tagged Jet + soft muon path: 

```
hltGetConfiguration --cff --offline /dev/CMSSW_11_2_0/GRun --paths HLTriggerFirstPath,HLTriggerFinalPath --unprescale > HLT_6AprV5_cff.py
hltGetConfiguration --cff /users/pdas/Tutorial2019/6Apr/AK8JetPath/V5 \
--globaltag auto:run3_mc_GRun   \
--unprescale >> HLT_6AprV5_cff.py
```
 You then need to modify HLT_6AprV5_cff.py : search for these two lines, they appear twice in the config file. Comment out their second appearance (not the first one). 
 ```
#/users/pdas/Tutorial2019/6Apr/AK8JetPath/V5 (CMSSW_11_2_0)   (already commented)                                                                                                                                       
#import FWCore.ParameterSet.Config as cms                                                                                                                                                                   
#fragment = cms.ProcessFragment( "HLT" )  
 ```
  
Copy the files to the python subdirectory
```
cp HLT_6AprV5_cff.py ../python/.
```

 Produce a HLT2_HLT.py config file to be run with cmsRun HLT2_HLT.py. We will need to rerun HLT using RAW data and and then access offline information from MINIAOD. In order to do that, we need to specify the MINIAOD file as the main input file and the list of all its parent RAW files as secondary input files: 
```
cmsDriver.py HLT2 --step=HLT:6AprV5 --era=Run3  --mc --conditions auto:run3_mc_GRun --filein file:/hdfs/store/user/pdas/Run3_ggHBB/MiniAOD/RunIIAutumn18MiniAOD_part0_75764563_75769863.root  --secondfilein file:/hdfs/store/user/pdas/Run3_ggHBB/DR/RunIIAutumn18DRPremix_step1_part0_75764563_75769863.root  --processName=HLT2 -n 10 --no_exec
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
For running L1boosted emulator, add:
```
process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')
process.load(‘L1Trigger.L1TCaloLayer1.uct2016EmulatorDigis_cfi’)

process.rawtodigi = cms.EndPath(process.l1tCaloLayer1Digis * process.uct2016EmulatorDigis)
```
Replace the line:
```
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])
```
by:
```
process.schedule.extend([process.endjob_step, process.rawtodigi, process.demo_step])
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
