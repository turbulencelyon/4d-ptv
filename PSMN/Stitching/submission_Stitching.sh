#! /bin/bash                                                                                                                                                                                                                                               

############  Parameters to tune  ######################
ManipName="Ra4.6e9_peudense2_6" #"Ra1.51e10_peudense_6"
dfmax=60              # Number of maximum tolerated missing frame 
dxmax=0.7                           # maximum distance between rays to consider a match
dvmax=0.35                               # minimum trajectory length
lmin=5                           # number of points used to predict next particle position
first=1                              # First frame of the experiment                                                                                           
last=36000                              # Last frame of the experiment
NbFramePerJobStitching=6000              # Number of frame per job for stitching. Has to be a multiple of NbFramePerJobTracking

Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where tge file track_x_x.h5 are 
Session_OUTPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the StitchA_x_x.h5 will be  

CompileFileDir="/Xnfs/convection/Lagrangien_David_Dumont/Toolbox/4d-ptv/PSMN/Stitching" #Directory where the file "run_Submission_Stitching.sh" is
LOG_path="$Session_OUTPUT/Processed_DATA/$ManipName/Parallel/Stitching/LOG"	#Log directory
OUT_path="$Session_OUTPUT/Processed_DATA/$ManipName/Parallel/Stitching/OUT"	#Matlab output directory
##############################################

export OMP_NUM_THREADS=1
for start in `seq $first $NbFramePerJobStitching $last`;
do
    stop=$(( start + NbFramePerJobStitching - 1 ));
    if [ $stop -gt $last ];
    then
        stop=$last;
    fi
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/stitching_$start-$stop.log -o $OUT_path/stitching_$start-$stop.log -q CLG5218deb192* -N Stitch-$start-$stop $CompileFileDir/run_submission_Stitching.sh $MCRROOT "$ManipName" "$start" "$stop" "$dfmax" "$dxmax" "$dvmax" "$lmin" "$Session_INPUT" "$Session_OUTPUT"

done;
