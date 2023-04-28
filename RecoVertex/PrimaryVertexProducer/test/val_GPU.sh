#!/bin/bash

eval `scramv1 runtime -sh`

cd ../../../
clear
cmsenv
#scram b clean
scram b -j 12
cd RecoVertex/PrimaryVertexProducer/test
cmsenv
cmsRun vertexTest.py

harvestTrackValidationPlots.py test_dqm_gpu.root -o gpu.root
rm -r plots/
makeTrackValidationPlots.py gpu.root cpu.root --png --extended
