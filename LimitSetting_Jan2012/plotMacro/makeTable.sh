#!/bin/bash



#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120121.beforePreapproval/singleChannel
#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120121.beforePreapproval/multiChannel

# for old and new comparison
#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120124.afterPreapproval/singleChannel

# check bin sensitivity with average backgrounds
#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120127/multiChannel

bino=$1
jet=$2
channel=$3

#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120131/${channel}
#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120203/${channel}
#dir=/uscms_data/d3/lpcpjm/jobs/20120203/${channel}
#dir=/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits.2012Jan/20120206/${channel}
dir=/uscms_data/d3/lpcpjm/jobs/20120209/${channel}

outfile=${bino}_${jet}.table

[ -e $outfile ] && rm $outfile
touch $outfile

prefix=$bino
[ $bino == "bino" ] && prefix="bino_mS"

for file in `dir -d $dir/${prefix}*${jet}*.result.txt`
do
    mS=`grep "mS :" $file | cut -d : -f 2`

    check=x$mS
    [ "$check" == "x" ] && continue

    mG=`grep "mG :" $file | cut -d : -f 2`
    mN=`grep "mN :" $file | cut -d : -f 2`
    acc=`grep "acc :" $file | cut -d : -f 2`
    xsec=`grep "xsecValue :" $file | cut -d : -f 2`
    xsecPDFError=`grep "xsecPDFError :" $file | cut -d : -f 2`
    xsecRSErrorNeg=`grep "xsecRSErrorNeg :" $file | cut -d : -f 2`
    xsecRSErrorPos=`grep "xsecRSErrorPos :" $file | cut -d : -f 2`

    obsLimit=`grep "CLs observed =" $file | cut -d = -f 2`
    expLimit=`grep "CLs expected =" $file | cut -d = -f 2`
    exp_m1s=`grep "CLs expected m1sigma =" $file | cut -d = -f 2`
    exp_m2s=`grep "CLs expected m2sigma =" $file | cut -d = -f 2`
    exp_p1s=`grep "CLs expected p1sigma =" $file | cut -d = -f 2`
    exp_p2s=`grep "CLs expected p2sigma =" $file | cut -d = -f 2`

    check=x$acc
    [ "$check" == "x 0" ] && continue

    check=x$exp_p2s
    [ "$check" == "x " ] && continue

    #echo "mS mG mN acc xsec xsecPDFError xsecRSErrorNeg xsecRSErrorPos obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s"

    echo "$mS $mG $mN $acc $xsec $xsecPDFError $xsecRSErrorNeg $xsecRSErrorPos $obsLimit $expLimit \
$exp_m1s $exp_m2s $exp_p1s $exp_p2s" >> $outfile

done

