#! /bin/bash                                                                                                                                                                                                                                               

############  Parameters to tune  ######################
ManipName="Ra1.51e10_peudense_6"
dfmax=60              # Number of frame per job for parallel matching
dxmax=0.5                           # maximum distance between rays to consider a match
dvmax=0.35                               # minimum trajectory length
lmin=5                           # number of points used to predict next particle position
first=35401                               # First frame of the experiment                                                                                           
last=36000                              # Last frame of the experiment
NbFramePerJobTracking=600              # Number of frame per job for tracking. Has to be a multiple of NbFramePerJob

Session_INPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where tge file track_x_x.h5 are 
Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the StitchA_x_x.h5 will be  

CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Stitching" #Directory where the file "run_Submission_Matlab_Stitching.sh" is
LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Stitching/LOG"	#Log directory
OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Stitching/OUT"	#Matlab output directory
##############################################

export OMP_NUM_THREADS=1
for start in `seq $first $NbFramePerJobTracking $last`;
do
    stop=$(( start + NbFramePerJobTracking - 1 ));
    if [ $stop -gt $last ];
    then
        stop=$last;
    fi
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/stitching_$start-$stop.log -o $OUT_path/stitching_$start-$stop.log -q piv_debian* -P PIV -N job-$start-$stop $CompileFileDir/run_Submission_Matlab_Stitching.sh $MCRROOT "$ManipName" "$start" "$stop" "$dfmax" "$dxmax" "$dvmax" "$lmin" "$Session_INPUT" "$Session_OUTPUT"

done;
