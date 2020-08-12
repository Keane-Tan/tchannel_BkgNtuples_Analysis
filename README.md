# tchannel_BkgNtuples_Analysis

The ntuples used in this analysis is stored in "root://cmseos.fnal.gov//store/user/keanet/CondorOutput/tchannel/BkgNtuples/MCRoot/". These ntuples were made from the nutples stored in "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV17/". The difference is that the ntuples used in this analysis have reduced amount of content in the tree, and each ntuple is a combination of multiple ntuples from the ones in Run2ProductionV17. For example, QCD16_Pt_80to120 in MCRoot is the combination of Summer16v3.QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8_i_RA2AnalysisTree.root, where i goes from 0 to 56 from Run2ProductionV17.

The output root files of this analysis scripts are stored in "root://cmseos.fnal.gov//store/user/keanet/CondorOutput/tchannelCom2020" (see RunCondor.sh). Each root file contains histograms of different kinematic distributions, which can be plotted for comparison between background and signal. 

To run the analysis script,
```
python FastSubmit.py
```

This should submit 16 jobs (4 backgrounds per period for 4 periods) to the lpc condor.
