#! /bin/bash

first=1
#last=2691
last=20
nperjob=2
export OMP_NUM_THREADS=1

for start in `seq $first $nperjob $last`;
do
    stop=$(( start + nperjob - 1 ));
    if [ $stop -gt $last ];
    then
        stop=$last;
    fi
    #if [ $start == 0 ];
    #then
    #    continue;
    #fi
    /usr/bin/qsub -cwd -V -e /home/eberna07/Test/job-$start-$stop.log -o /home/eberna07/Test/job-$start-$stop.log -q piv_debian -P PIV -N dr-$start-$stop /home/eberna07/Stage_EB_2020/4d-ptv/Matching/STMPython/stm.py "/home/eberna07/Stage_EB_2020/4d-ptv/Documentation/TestData/Processed_DATA/MyExperiment/rays2.dat" 1 20 2 0.2 400 400 250 2;
done;
