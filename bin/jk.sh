#!/bin/bash

display_usage() { 
    echo "Shell script for partitioning 2D and 3D data points into N^2 and N^3 regions."
    echo "Partitioning is performed by equal number of points (+/-1)."
    echo "Infile format 3-4 columns eg (x y z weight) or (RA DEC weight)"
    echo "Output format, same as input with extra last column of JK region."
    echo ""
    echo -e "\nUsage:\n jk.sh galfile ranfile N D \n" 
} 
# check whether user had supplied -h or --help . If yes display usage 
if [[ ( $# == "--help") ||  $# == "-h" ]] 
    then 
    display_usage
    exit 0
    fi 

if [  $# -le 1 ] 
then 
    display_usage
    exit 1
fi 


x=$1
y=$2
echo "working on " $x $y
num_files=$3
D=$4

if [  $D == 3 ] 
then 
    echo "Working with 3D data" 
    echo "chopping into" $((num_files*num_files*num_files)) "regions"
elif [ $D == 2 ]
then
    echo "Working with 2D data" 
    echo "chopping into" $((num_files*num_files)) "regions"
else
    echo "Unknown Value of D. Either 2 or 3"
    exit 1
fi


awk '{print $0, 1, NR}' $x > $x.gal
awk '{print $0, 0, NR}' $y > $y.ran
cat $x.gal $y.ran > new.txt
sort -gk1 new.txt > new.txt.sort

total_lines=$(wc -l < new.txt.sort)
echo number of lines: $total_lines

((lines_per_file = (total_lines + num_files - 1) / num_files))
echo lines per file: $lines_per_file

#######################################################################
split -d -l ${lines_per_file} new.txt.sort firstsplit.
echo performing x-axis split....
n=0
while [ $n -le $((num_files - 1)) ]
do
awk -v lin="$n" '{print $0, lin }' firstsplit.0$n > firstsplit2.$n
n=$((n + 1))
done

#######################################################################
echo performing y-axis split....
n=0
while [ $n -le $((num_files - 1)) ]
do
sort -gk2 firstsplit2.$n > firstsplit2.${n}.sort
total_lines=$(wc -l < firstsplit2.${n}.sort)

((lines_per_file = (total_lines + num_files - 1) / num_files))

split -d -l ${lines_per_file} firstsplit2.${n}.sort secondsplit.${n}.
m=0
while [ $m -le $((num_files - 1)) ]
do
awk -v lin="$m" '{print $0, lin }' secondsplit.$n.0$m > secondsplit2.$n.$m
m=$((m + 1))
done
n=$((n + 1))
done
cat secondsplit2* > secondsplit

if [  $D == 3 ] 
then 
#######################################################################
echo performing z-axis split....
n=0
while [ $n -le $((num_files - 1)) ]
do
m=0
while [ $m -le $((num_files - 1)) ]
do
sort -gk3 secondsplit2.${n}.${m} > secondsplit2.${n}.${m}.sort
total_lines=$(wc -l < secondsplit2.${n}.${m}.sort)

((lines_per_file = (total_lines + num_files - 1) / num_files))

split -d -l ${lines_per_file} secondsplit2.${n}.${m}.sort thirdsplit.${n}.${m}.
p=0
while [ $p -le $((num_files - 1)) ]
do
awk -v lin="$p" '{print $0, lin }' thirdsplit.$n.$m.0$p > thirdsplit2.$n.$m.$p
p=$((p + 1))
done
m=$((m + 1))
done
n=$((n + 1))
done

cat thirdsplit2.* > thirdsplit
fi

#######################################################################



if [  $D == 3 ] 
then 
    total_lines=$(wc -l < thirdsplit)
    echo number of lines: $total_lines
    awk -v nn=$num_files '($5==1){print $1,$2,$3,$4,(1 + $7 + nn*$8 + nn*nn*$9)}' thirdsplit > ${x}.jk
    awk -v nn=$num_files '($5==0){print $1,$2,$3,$4,(1 + $7 + nn*$8 + nn*nn*$9)}' thirdsplit > ${y}.jk
else
    total_lines=$(wc -l < secondsplit)
    echo number of lines: $total_lines
    awk -v nn=$num_files '($3==1){print $1,$2,(1 + $5 + nn*$6)}' secondsplit > ${x}.jk
    awk -v nn=$num_files '($3==0){print $1,$2,(1 + $5 + nn*$6)}' secondsplit > ${y}.jk 
fi

rm -fr firstsplit* secondsplit* thirdsplit* new.txt new.txt.sort $x.gal $y.ran
