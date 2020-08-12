#!/bin/bash

mkdir kAnalysis
cd kAnalysis

source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a tcsh script, use .csh instead of .sh
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 project CMSSW CMSSW_10_2_20`
cd CMSSW_10_2_20/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
echo "CMSSW: "$CMSSW_BASE

OUTDIR=root://cmseos.fnal.gov//store/user/keanet/CondorOutput/tchannelCom2020

# transfering macro files
xrdcp -s $OUTDIR/tchannelCom.tgz .
tar -xf tchannelCom.tgz # untarring the tar ball
echo "After cd-ing into CMSSW_10_2_20/src and ls-ing"
ls
cd tchannelCom
echo "After cd-ing into N1Analysis and ls-ing"
ls
echo "python tchannelCom.py ${1} ${2}"
python tchannelCom.py ${1} ${2}

echo "xrdcp -f output files"
xrdcp -f tCom_${2}${1}.root $OUTDIR/tCom_${2}${1}.root

# Remove everything
cd ${_CONDOR_SCRATCH_DIR}
rm -rf kAnalysis
