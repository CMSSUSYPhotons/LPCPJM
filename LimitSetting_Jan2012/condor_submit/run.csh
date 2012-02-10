#!/bin/tcsh

set sec=$1
set bino=$2

set WORK_DIR = `pwd`
set SRC_DIR = /uscms/home/dwjang/rel/428/src/ra3
#set DATA_DIR = /uscms/home/dwjang/work/jobs/datacards/20120124_combinedBG/singleChannel
#set DATA_DIR = /uscms/home/dwjang/work/jobs/datacards/20120124/multiChannel
set DATA_DIR = /uscms/home/dwjang/work/jobs/datacards/20120131/multiChannel
set TEMP_DIR = ${WORK_DIR}/temp

cd $SRC_DIR
source /uscmst1/prod/sw/cms/cshrc prod
eval `scramv1 runtime -csh`

cd $WORK_DIR

mkdir -p $TEMP_DIR

cd $TEMP_DIR
cp $SRC_DIR/condor_submit/limit .
ln -s $SRC_DIR/xsecdat .

if ($bino == "bino" || $bino == "wino") then
    set iter=0
    set iS = 400
    while ($iS <= 2000)
	set iG = 400
	while ($iG <= 2000)
	    if($iter == $sec) then
		set mS = $iS
		set mG = $iG
	    endif
	    @ iter++
	    @ iG = $iG + 80
	end
	@ iS = $iS + 80
    end
    set label=${bino}_mS${mS}_mG${mG}_mN375
endif

if ($bino == "bino_mNScan") then
    set iter=0
    set iN = 150
    while ($iN <= 1050)
	set iG = 160
	while ($iG <= 2000)
	    if($iter == $sec) then
		set mN = $iN
		set mG = $iG
	    endif
	    @ iter++
	    @ iG = $iG + 80
	end
	@ iN = $iN + 100
    end
    set label=${bino}_mS2500_mG${mG}_mN${mN}
endif

set logfile = $WORK_DIR/log.${bino}_sec${sec}_${label}
touch $logfile


echo                                                     >>& $logfile
echo "-------------------------------------------------" >>& $logfile
echo "section $sec nojet start..."                       >>& $logfile
echo "-------------------------------------------------" >>& $logfile
echo                                                     >>& $logfile
mkdir $TEMP_DIR/nojet
cd $TEMP_DIR/nojet
$TEMP_DIR/limit $DATA_DIR/${label}_nojet.dat   >>& $logfile
mv *.result.txt $WORK_DIR

echo                                                     >>& $logfile
echo                                                     >>& $logfile
echo "-------------------------------------------------" >>& $logfile
echo "section $sec 1jet start..."                        >>& $logfile
echo "-------------------------------------------------" >>& $logfile
echo                                                     >>& $logfile
mkdir $TEMP_DIR/1jet
cd $TEMP_DIR/1jet
$TEMP_DIR/limit $DATA_DIR/${label}_1jet.dat    >>& $logfile
mv *.result.txt $WORK_DIR

echo                                                     >>& $logfile
echo                                                     >>& $logfile
echo "-------------------------------------------------" >>& $logfile
echo "moving result files and cleaning..."               >>& $logfile
echo "-------------------------------------------------" >>& $logfile

cd $WORK_DIR

rm -rf $TEMP_DIR

