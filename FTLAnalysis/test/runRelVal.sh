#!/bin/sh

runSample()
{
    python drawTrackPerformance.py --inputDir=$1  --pattern='DumpHits*' --output=$1/drawTrackPerformance.root --layout=barzflat 
}

runSample_noMTD()
{
    python drawTrackPerformance.py --inputDir=$1  --pattern='DumpHits_noMTD*' --output=$1/drawTrackPerformance.root --layout=barzflat 
}

### V8
#runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10_pythia8/runHitsRelValSingleMuFlatPt0p7to10pythia8_v8/190112_165725/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/runHitsRelValSinglePiFlatPt0p7to10pythia8cfi_v8/190112_165643/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValDYToLL_M_50_14TeV/runHitsRelValDYToLLM5014TeV_v8/190112_170107/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10/runHitsRelValSingleMuFlatPt0p7to10_v8/190112_165809/0000
#runSample /eos/cms/store/user/meridian/MTD/DYToLL_M-50_14TeV_pythia8/runHitsDYToLLM-5014TeVpythia8_v8/190112_170155/0000

### V7
runSample /eos/cms/store/user/meridian/MTD/DYToLL_M-50_14TeV_pythia8/runHitsDYToLLM-5014TeVpythia8_v7/190112_192028/0000
runSample /eos/cms/store/user/meridian/MTD/RelValDYToLL_M_50_14TeV/runHitsRelValDYToLLM5014TeV_v7/190112_191946/0000
runSample /eos/cms/store/user/meridian/MTD/RelValMinBias_14TeV/runHitsRelValMinBias14TeV_v7/190112_191903/0000
runSample /eos/cms/store/user/meridian/MTD/RelValSingleKaonFlatPt_0p7to10/runHitsRelValSingleKaonFlatPt0p7to10_v7/190112_191738/0000
runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10/runHitsRelValSingleMuFlatPt0p7to10_v7/190112_191655/0000
runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10_pythia8/runHitsRelValSingleMuFlatPt0p7to10pythia8_v7/190112_191612/0000
runSample /eos/cms/store/user/meridian/MTD/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/runHitsRelValSinglePiFlatPt0p7to10pythia8cfi_v7/190112_191423/0000
runSample /eos/cms/store/user/meridian/MTD/RelValSingleProtonFlatPt_0p7to10/runHitsRelValSingleProtonFlatPt0p7to10_v7/190112_191304/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/DYToLL_M-50_14TeV_pythia8/runHitsDYToLLM-5014TeVpythia8_v7_noMTD/190113_070207/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValDYToLL_M_50_14TeV/runHitsRelValDYToLLM5014TeV_v7_noMTD/190113_070123/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValMinBias_14TeV/runHitsRelValMinBias14TeV_v7_noMTD/190113_070038/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleKaonFlatPt_0p7to10/runHitsRelValSingleKaonFlatPt0p7to10_v7_noMTD/190113_065908/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10/runHitsRelValSingleMuFlatPt0p7to10_v7_noMTD/190113_065824/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10_pythia8/runHitsRelValSingleMuFlatPt0p7to10pythia8_v7_noMTD/190113_065739/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/runHitsRelValSinglePiFlatPt0p7to10pythia8cfi_v7_noMTD/190113_065654/0000
runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleProtonFlatPt_0p7to10/runHitsRelValSingleProtonFlatPt0p7to10_v7_noMTD/190113_065610/0000

### V9

#runSample /eos/cms/store/user/meridian/MTD/DYToLL_M-50_14TeV_pythia8/runHitsDYToLLM-5014TeVpythia8_v9/190113_173727/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValDYToLL_M_50_14TeV/runHitsRelValDYToLLM5014TeV_v9/190113_173644/0000
##runSample /eos/cms/store/user/meridian/MTD/RelValMinBias_14TeV/runHitsRelValMinBias14TeV_v9/190113_173602/0000
##runSample /eos/cms/store/user/meridian/MTD/RelValSingleKaonFlatPt_0p7to10/runHitsRelValSingleKaonFlatPt0p7to10_v9/190113_173436/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10/runHitsRelValSingleMuFlatPt0p7to10_v9/190113_173353/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10_pythia8/runHitsRelValSingleMuFlatPt0p7to10pythia8_v9/190113_173311/0000
#runSample /eos/cms/store/user/meridian/MTD/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/runHitsRelValSinglePiFlatPt0p7to10pythia8cfi_v9/190113_173228/0000
##runSample /eos/cms/store/user/meridian/MTD/RelValSingleProtonFlatPt_0p7to10/runHitsRelValSingleProtonFlatPt0p7to10_v9/190113_173145/0000
#
#runSample_noMTD /eos/cms/store/user/meridian/MTD/DYToLL_M-50_14TeV_pythia8/runHitsDYToLLM-5014TeVpythia8_v9_noMTD/190113_172852/0000
#runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValDYToLL_M_50_14TeV/runHitsRelValDYToLLM5014TeV_v9_noMTD/190113_172810/0000
##runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValMinBias_14TeV/runHitsRelValMinBias14TeV_v9_noMTD/190113_172728/0000
##runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleKaonFlatPt_0p7to10/runHitsRelValSingleKaonFlatPt0p7to10_v9_noMTD/190113_172603/0000
#runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10/runHitsRelValSingleMuFlatPt0p7to10_v9_noMTD/190113_172519/0000
#runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10_pythia8/runHitsRelValSingleMuFlatPt0p7to10pythia8_v9_noMTD/190113_172437/0000
#runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/runHitsRelValSinglePiFlatPt0p7to10pythia8cfi_v9_noMTD/190113_172354/0000
##runSample_noMTD /eos/cms/store/user/meridian/MTD/RelValSingleProtonFlatPt_0p7to10/runHitsRelValSingleProtonFlatPt0p7to10_v9_noMTD/190113_172311/0000

### V6
#runSample /eos/cms/store/user/meridian/MTD/RelValSingleMuFlatPt_0p7to10_pythia8/runHitsRelValSingleMuFlatPt0p7to10pythia8_v6/190112_080128/0000/

### V5
#runSample /eos/cms/store/user/meridian/MTD/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/runHitsRelValSinglePiFlatPt0p7to10pythia8cfi_v5/190111_145314/0000/