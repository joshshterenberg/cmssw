#!/bin/bash

eval `scramv1 runtime -sh`

cd ../../../
clear
cmsenv
#scram b clean
scram b -j 12 USER_CXXFLAGS="-g"
cd RecoVertex/PrimaryVertexProducer/test
cmsenv
cmsRun vertexTest.py
