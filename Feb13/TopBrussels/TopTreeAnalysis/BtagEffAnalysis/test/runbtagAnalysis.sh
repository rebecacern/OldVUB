#!/bin/bash   

echo "BASH::STARTING BASH SCRIPT"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

dobin="ninbin" # not used
dobin="noutbin" # not used

for chicut in 6; do # GOOD

#for chicut in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do

#for chicut in 15 16 17 18 19; do
#
iFbias=-1 # only MC

### For testing Fbias=true
#for iFbias in 0 1 2 3 4 5 6 7 8 9 10; do
#chicut=23
####################    

    #sample=0
    #sample=1
    #sample=2
    #sample=3
    #sample=4
    #sample=5
    #sample=6

    OUTBIN=0 # not used
    INBIN=0 # not used
    DIRNAME="a" # not used


    #chiSqCut=18
    chiSqCut=$chicut
    echo chiSqCut = $chiSqCut

    DIRNAME="TTrees/DATA_$chiSqCut"
    ####################    
   
    if [ $dobin = "inbin" ]; then
	echo "BASH::DOING INBIN"
	#INBIN=0
	INBIN=$bin
	echo inbin = $INBIN
	#DIRNAME=../outRootFiles_29/TT_MG_ST_W_Z_1fb_JES+10%_chisq_18_inBin_$INBIN
	#DIRNAME=../outRootFiles_29/TT_MG_ST_W_Z_1fb_newF_JES+0%_ptbin_chisq_18_inBin_$INBIN
	#DIRNAME=../outRootFiles_29/TTAll_ptbinmore_chisq_18_inBin_$INBIN
	#DIRNAME=../outRootFiles_31/TT_MG_ST_W_Z_ptbinmore_chisq_18_inBin_$INBIN
	
	OUTBIN=0
    fi    

    if [ $dobin = "outbin" ]; then
	OUTBIN=$bin
	echo outbin = $OUTBIN
	#DIRNAME=../outRootFiles_29/TT_MG_ST_W_Z_1fb_chisq_18_bin3_$OUTBIN
	echo dirname = $DIRNAME
	INBIN=0
    fi

    mkdir $DIRNAME
    
    cp runbtagAnalysis.sh $DIRNAME
    cp mybTagMeasurement.C $DIRNAME
    
    
    #FILENAME="info_bias_"$i"_and_bg_"$j".txt" 
    FILENAME="info$(date +%s).txt" 
    
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../..:../../TMVA/lib/:/Users/michael/lib:/user/mmaes/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:../..:../../TMVA/lib/:/Users/michael/lib:/user/mmaes/lib

    #echo g++ -m64 -L `pwd` -l Toto -l bTagAnalysis -l TopTreeAna -I `root-config --incdir` `root-config --libs` mybTagMeasurement.C -o mybTagMeasurement 
    #echo g++ -m64 -L../../ -L../../TMVA/lib/ -lTopTreeAnaContent -lTopTreeAna -I `root-config --incdir` `root-config --libs` mybTagMeasurement.C -o mybTagMeasurement 
    #echo g++ -m64 -L../../ -lBtagAnalysis -I `root-config --incdir` `root-config --libs` mybTagMeasurement.C -o mybTagMeasurement 

    g++ -m64 -I../../../ -I../../ -I.. -L/user/mmaes/lib -L/Users/michael/lib -lBtagAnalysis53 -I `root-config --incdir` `root-config --libs` mybTagMeasurement.C -o mybTagMeasurement
    #The value on the 4th position will tell you how much of the background you will remove (exact value from a list in mybTagMeasurement)

    ### For testing Fbias=true
	
	if [ $? == "0" ];then	
		echo ./mybTagMeasurement 2 $chiSqCut $iFbias -1 -1 $INBIN $OUTBIN $* | tee $DIRNAME/$FILENAME 
    ./mybTagMeasurement 2 $chiSqCut 0 $iFbias -1 -1 $INBIN $OUTBIN $* | tee $DIRNAME/$FILENAME 


fi
	####################    

    ### For runnin individual samples
    #echo ./mybTagMeasurement 2 $chiSqCut -1 -1 $sample $INBIN $OUTBIN | tee $DIRNAME/$FILENAME 
    #./mybTagMeasurement 2 $chiSqCut -1 -1 $sample $INBIN $OUTBIN | tee $DIRNAME/$FILENAME 
    ####################    


    #echo ./mybTagMeasurement 2 $chiSqCut -1 -1 -1 $INBIN $OUTBIN | tee $DIRNAME/$FILENAME 
    #./mybTagMeasurement 2 $chiSqCut -1 -1 -1 $INBIN $OUTBIN | tee $DIRNAME/$FILENAME 

    #echo ./mybTagMeasurement 2 $chiSqCut -1 -1 0 $INBIN $OUTBIN | tee $DIRNAME/$FILENAME 
    #./mybTagMeasurement 2 $chiSqCut -1 -1 0 $INBIN $OUTBIN | tee $DIRNAME/$FILENAME 
    
done

echo "BASH::FINISHING BASH SCRIPT"


