#!/bin/bash
#
# run normal mode program to find normal mode solution tables
#

if [ $# != 2 ]; then
echo "### Usage of $0 ###"
echo "% run.sh 'input file for solid_rk' 'output dispersion table file'"
echo "### End of Usage of $0 ###"
exit
fi
#
SOLID_RK_IN=$1
base=`basename $1`
indir=`dirname $1`
outdir=`dirname $2`
SOLID_RK_OUT=$outdir/solid_rk.out${base#solid_rk.in}
MODE_TABLE_FILE=$2

if [ ! -f  $SOLID_RK_IN ]; then
echo "input file [ $SOLID_RK_IN ] does not exist."
exit
fi

if [ -f $MODE_TABLE_FILE ]; then
echo "output file [ $MODE_TABLE_FILE ] already exists." 
exit
fi

#
programs/solid_rk < $SOLID_RK_IN > $SOLID_RK_OUT
#
# sort the normal mode solutions ( agular order, eigenfrequency, period, phase velocity, group velocity )
# and make a table
# 
echo "#l omg period V_phase V_group PREM 4km ocean" > $MODE_TABLE_FILE
paste -d" " <(grep found $SOLID_RK_OUT | awk 'BEGIN{FS="="}{print $2,$3,$4}' | awk '{print $1,$3,$5}' )\
 <( grep Vc $SOLID_RK_OUT | awk 'BEGIN{FS="="}{print $3,$4}' | awk '{print $1,$3}') \
>> $MODE_TABLE_FILE
