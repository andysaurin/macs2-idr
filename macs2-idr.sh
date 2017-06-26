#!/bin/bash
#
# This script calls IDR peaks on 2 ChIP replicates using Input as control
# It is advised to check the macs_output directory for MACS2 problems
#
# Author: Andy Saurin
# Email: andrew.saurin at univ-amu.fr
#
#
###########################################################################
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	Copyright 2014 Andy Saurin
#
###########################################################################

INSTALL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

BATCH_CONSISTENCY_ANALYSIS=${INSTALL_DIR}"/idrCode/batch-consistency-analysis.r"

if ! [[ -f $BATCH_CONSISTENCY_ANALYSIS ]]; then
	echo "The batch-consistency-analysis.r Rscript cannot be found."
	echo "Please edit line 31 of ( ${BASH_SOURCE[0]} ) to point to where the batch-consistency-analysis.r can be found."
	echo ""
	echo "You can download batch-consistency-analysis.r at https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz?attredirects=0"
	echo "Exiting."
	exit 1
fi

COMPARE_IDR_PEAKS=${INSTALL_DIR}"/compare_idr_peaks.R"


command -v macs2 >/dev/null 2>&1 || { echo >&2 "$0 requires MACS2, but it's not installed or in your path.  Aborting."; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo >&2 "$0 requires Rscript, but it's not installed or in your path.  Aborting."; exit 1; }

PWD=`pwd`

usage()
{
cat << EOF
usage: $0 options

This script will call IDR peaks on two ChIP replicats using MACS2, with 1 control input.

OPTIONS:
   -f      Factor name (Required. This will be used for output file naming)
   -1     The ChIP Rep 1 BAM file (Required.)
   -2     The ChIP Rep 2 BAM file (Required.)
   -c      The Control (Input) BAM file (Required.)
   -o      Output directory (default: ./peak_calling)
   -p      Peak calling type. One of:
               narrow (Default) - macs2 narrowPeak calling
               broad - macs2 broadPeak calling
   -g      Genome to use. One of:
               dm - Drosophila (Default)
               mm - Mouse
               ce - C. elegans
               hs - Human
   -m      Optional commands to pass through to MACS2. Provide as a QUOTED option eg: -m "--bw 400 --mfold 4 100"
   -h      Show this message

EOF
}

FACTOR=
CHIP1=
CHIP2=
INPUT=
OUTDIR=`pwd`/peak_calling
MACS_OPTIONS=
PEAKTYPE='narrow'
GENOME='dm'

while getopts "hf:1:2:c:o:p:g:m:" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FACTOR=$OPTARG
             ;;
         1)
             CHIP1=`readlink -f $OPTARG`
             ;;
         2)
             CHIP2=`readlink -f $OPTARG`
             ;;
         c)
             INPUT=`readlink -f $OPTARG`
             ;;
         o)
             OUTDIR=`readlink -f $OPTARG`
             ;;
         p)
             PEAKTYPE=$OPTARG
             ;;
         g)
             GENOME=$OPTARG
             ;;
         m)
             MACS_OPTIONS=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


# Remove quotes from macs options
if ! [[ -z $MACS_OPTIONS ]]; then
	MACS_OPTIONS="${MACS_OPTIONS%\"}"
	MACS_OPTIONS="${MACS_OPTIONS#\"}"
	MACS_OPTIONS="${MACS_OPTIONS%\'}"
	MACS_OPTIONS="${MACS_OPTIONS#\'}"
fi



if [[ -z $FACTOR ]] || [[ -z $CHIP1 ]] || [[ -z $CHIP2 ]] || [[ -z $INPUT ]]; then
     usage
     exit 1
fi


if ! [[ -f $CHIP1 ]]; then
	usage

	echo "No such BAM file for ChIP Rep1: "${CHIP1}
	exit 1
fi

if ! [[ -f $CHIP2 ]]; then
	usage

	echo "No such BAM file for ChIP Rep2: "${CHIP2}
	exit 1
fi

if ! [[ -f $INPUT ]]; then
	usage

	echo "No such BAM file for Input: "${INPUT}
	exit 1
fi

if ! [[ -d $OUTDIR ]]; then
	mkdir -p $OUTDIR || { echo "Cannot create output directory $OUTDIR"; exit 1; }
fi

if  [[ $GENOME != 'dm' ]] &&  [[ $GENOME != 'mm' ]] &&  [[ $GENOME != 'ce' ]] &&  [[ $GENOME != 'hs' ]]; then
	usage

	echo "-g can be one of either dm, mm, ce, hs. Default is 'dm'"
	exit 1
fi

if  [[ $PEAKTYPE != 'narrow' ]] &&  [[ $PEAKTYPE != 'broad' ]]; then
	usage

	echo "-p is either 'narrow' or 'broad'. -p narrow is default"
	exit 1
fi

if  [[ $PEAKTYPE == 'narrow' ]]; then
	MACSPEAKTYPE=''
	CALLSUMMITS='--call-summits'
fi
if  [[ $PEAKTYPE == 'broad' ]]; then
	MACSPEAKTYPE='--broad'
	CALLSUMMITS=''
fi


rotateCursor()
{
  case $toggle
  in
    1)
      echo -n $1" \ "
      echo -ne "\r"
      toggle="2"
    ;;

    2)
      echo -n $1" | "
      echo -ne "\r"
      toggle="3"
    ;;

    3)
      echo -n $1" / "
      echo -ne "\r"
      toggle="4"
    ;;

    *)
      echo -n $1" - "
      echo -ne "\r"
      toggle="1"
    ;;
  esac
}




mkdir -p ${OUTDIR}/${FACTOR}/idr
mkdir -p ${OUTDIR}/${FACTOR}/macs_output

MACS_OUTPUT=${OUTDIR}/${FACTOR}/macs_output
cd ${OUTDIR}/${FACTOR}

MACS_LOG_REP1=${MACS_OUTPUT}/log_rep1.txt
MACS_LOG_REP2=${MACS_OUTPUT}/log_rep2.txt


# MACS2 on Rep1
echo "Calling peaks on ChIP 1"

macs2 callpeak -t $CHIP1 -c $INPUT \
  $MACS_OPTIONS \
  --outdir=$MACS_OUTPUT $CALLSUMMITS $MACSPEAKTYPE -p 0.1 -g $GENOME \
  -n ${FACTOR}_Rep1_p0.1 > $MACS_LOG_REP1 2>&1 &

pid=$!
trap "kill $pid 2> /dev/null" EXIT
while kill -0 $pid 2> /dev/null; do
	rotateCursor
	sleep 1
done
# Disable the trap on a normal exit.
trap - EXIT
echo -n "Done. "

NUM_WARNINGS_REP1=`grep 'WARNING ' ${MACS_LOG_REP1} | wc -l`
if (( $NUM_WARNINGS_REP1 > 0 )); then
	echo "The following MACS2 warnings were issued:"
	echo `grep 'WARNING ' ${MACS_LOG_REP1}`
else
	echo "No problems were encountered :)"
fi

echo ""

# MACS2 on Rep2
echo "Calling peaks on ChIP 2"

macs2 callpeak -t $CHIP2 -c $INPUT \
  $MACS_OPTIONS \
  --outdir=$MACS_OUTPUT $CALLSUMMITS $MACSPEAKTYPE -p 0.1 -g $GENOME \
  -n ${FACTOR}_Rep2_p0.1 > $MACS_LOG_REP2 2>&1 &

pid=$!
trap "kill $pid 2> /dev/null" EXIT
while kill -0 $pid 2> /dev/null; do
	rotateCursor
	sleep 1
done
# Disable the trap on a normal exit.
trap - EXIT
echo -n "Done. "

NUM_WARNINGS_REP2=`grep 'WARNING ' ${MACS_LOG_REP2} | wc -l`
if (( $NUM_WARNINGS_REP2 > 0 )); then
	echo "The following MACS2 warnings were issued:"
	echo `grep 'WARNING ' ${MACS_LOG_REP2}`
else
	echo "No problems were encountered :)"
fi

echo ""

if [[ $PEAKTYPE == 'broad' ]]; then
	FILETYPE='broadPeak';
	ISBROADPEAK=T
fi
if [[ $PEAKTYPE == 'narrow' ]]; then
	FILETYPE='narrowPeak';
	ISBROADPEAK=F
fi

if ! [[ -f ${OUTDIR}/${FACTOR}/macs_output/${FACTOR}_Rep1_p0.1_peaks.${FILETYPE} ]]; then
	echo "MACS2 failed on Rep1"
	exit 1
fi
if ! [[ -f ${OUTDIR}/${FACTOR}/macs_output/${FACTOR}_Rep2_p0.1_peaks.${FILETYPE} ]]; then
	echo "MACS2 failed on Rep2"
	exit 1
fi

echo "Running batch-consistency-analysis.r"


cd $PWD

Rscript $BATCH_CONSISTENCY_ANALYSIS \
  ${OUTDIR}/${FACTOR}/macs_output/${FACTOR}_Rep1_p0.1_peaks.${FILETYPE} \
  ${OUTDIR}/${FACTOR}/macs_output/${FACTOR}_Rep2_p0.1_peaks.${FILETYPE} \
  -1 ${OUTDIR}/${FACTOR}/idr/${FACTOR}_Rep1_inter_Rep2.MACS2 0 ${ISBROADPEAK} p.value &

pid=$!
trap "kill $pid 2> /dev/null" EXIT
while kill -0 $pid 2> /dev/null; do
	rotateCursor
	sleep 1
done
# Disable the trap on a normal exit.
trap - EXIT
echo "Done"

if ! [[ -f ${OUTDIR}/${FACTOR}/idr/${FACTOR}_Rep1_inter_Rep2.MACS2-overlapped-peaks.txt ]]; then
	echo "Batch consistency analysis failed"
	exit 1
fi

echo ""
echo -n "Creating IDR BED files"

for ((i=1;i<=10;i++)); do
   if [ $i -eq 10 ]
   then
      frac="0.1"
   else
      frac="0.0${i}"
   fi
   sed 's/\"//g' ${OUTDIR}/${FACTOR}/idr/${FACTOR}_Rep1_inter_Rep2.MACS2-overlapped-peaks.txt | \
    awk -F' ' '{if ( $11 <= '$frac' ) { print $0 } }' - | \
    awk 'BEGIN{x=1; frac=$frac}; { \
     if ($2 ~ /^chr/) { print $2"\t"$3"\t"$4"\tMACS2_IDR_peak_"x"\t"$5; ((x++)) } \
    }' - > ${OUTDIR}/${FACTOR}/idr/${FACTOR}.MACS2_IDR-${frac}_peaks.bed
    awk -v frac=$frac -v factor=$FACTOR 'BEGIN{x=1}; {  print $1"\t"$2"\t"$3"\t"factor"_MACS2_IDR-"frac"_peak_"x"\t"$5; ((x++)) }' ${OUTDIR}/${FACTOR}/idr/${FACTOR}.MACS2_IDR-${frac}_peaks.bed > ${OUTDIR}/${FACTOR}/${FACTOR}.MACS2_IDR-${frac}_peaks.bed && rm -f ${OUTDIR}/${FACTOR}/idr/${FACTOR}.MACS2_IDR-${frac}_peaks.bed
done


echo ""

echo "###################################"
echo "             SUMMARY              "
echo "                                  "

for ((i=1;i<=10;i++)); do
   if [ $i -eq 10 ]
   then
      frac="0.1"
   else
      frac="0.0${i}"
   fi
   if [[ -f ${OUTDIR}/${FACTOR}/${FACTOR}.MACS2_IDR-${frac}_peaks.bed ]]; then
      echo " IDR ${i}%: "`wc -l ${OUTDIR}/${FACTOR}/${FACTOR}.MACS2_IDR-${frac}_peaks.bed | awk '{print $1}' -`" peaks"
   fi
done

echo "                                  "
echo "                                  "
echo "###################################"

echo ""
echo "Output IDR BED files in ${OUTDIR}/${FACTOR}"

echo ""

echo "Generating graphs of different IDR peaks called"
echo ""

${COMPARE_IDR_PEAKS} -F ${FACTOR} -D ${OUTDIR}
