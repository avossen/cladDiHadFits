#!/bin/bash
rm doFits
CFLAGS=" -Wall -ggdb `root-config --cflags --libs`  -lMinuit -lRooFitCore -lRooFit"
echo cflags: $CFLAGS
#c++ $CFLAGS doAnalysis.cc -o doAnalysis
c++ $CFLAGS doFits.cc  -o doFits
#c++ $CFLAGS test2Sys.cc binning.cc sysDataGen.cc -o test2Sys
#c++ $CFLAGS doAnalysis2.cc -o doAnalysis2
#c++ $CFLAGS doSoBAnalysis.cc -o doSoBAnalysis
#c++ $CFLAGS doSoBAnalysis_upperCut.cc -o doSoBAnalysis_upperCut
#c++ $CFLAGS doSoBAnalysis_lowerCut.cc -o doSoBAnalysis_lowerCut
#c++ $CFLAGS doSoBAnalysisPionMomCut.cc -o doSoBAnalysisPionMomCut
#c++ $CFLAGS doSoBAnalysisPion1MomCut.cc -o doSoBAnalysisPion1MomCut
#c++ $CFLAGS doSoBAnalysisPion1LowerMomCut.cc -o doSoBAnalysisPion1LowerMomCut
#c++ $CFLAGS doSoBAnalysisPion2LowerMomCut.cc -o doSoBAnalysisPion2LowerMomCut
###c++ $CFLAGS doSoBAnalysisLogProbCutHuschle.cc binning.cc -o doSoBAnalysisLogProbCutHuschle
#c++ $CFLAGS doSoBAnalysisMDnCut.cc -o doSoBAnalysisDnCut
#c++ $CFLAGS doSoBAnalysisMinMDnCut.cc -o doSoBAnalysisMinDnCut
#c++ $CFLAGS doSoBAnalysisDeltaECut.cc binning.cc -o doSoBAnalysisDeltaECut
###c++ $CFLAGS doSoBAnalysisDeltaECutHuschle.cc binning.cc -o doSoBAnalysisDeltaECutHuschle
###c++ $CFLAGS doSoBAnalysisDeltaMCutHuschle.cc binning.cc -o doSoBAnalysisDeltaMCutHuschle
#c++ $CFLAGS doSoBAnalysis_DecayChannels.cc -o doSoBAnalysis_DecayChannels


#c++ $CFLAGS doOverlapCount.cc -o doOverlapCount

