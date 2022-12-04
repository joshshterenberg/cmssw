#!/bin/bash

clear
eval `scramv1 runtime -sh`

cmsenv
cd ../../../
scram b -j 12
cd RecoVertex/PrimaryVertexProducer/test

cmsRun vertexTest.py gpu=True # make sure this is what u want loser

harvestTrackValidationPlots.py test_dqm_gpu.root -o gpu.root

rm -r plots/
rm -r ~/plots/

makeTrackValidationPlots.py cpu_gauss.root gpu.root --png --extended
mv plots/ /eos/user/j/jshteren/www/plotsyeet2
