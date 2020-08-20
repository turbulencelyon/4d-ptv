#! /bin/bash                                                                                                                                                                                                                                               

############  Parameter to tune  ######################
NbFramePerJobMatching=20             # Number of frame per job for parallel matching
maxdist=0.4                         # maximum distance between rays to consider a match
lmin=5                                # minimum trajectory length
npriormax=5                           # number of points used to predict next particle position
manipname="Ra1.51e10_peudense_6"
first=401                               # First frame of the experiment                                                                                           
last=36000                           # Last frame of the experiment
NbFramePerJobTracking=5000             # Number of frame per job for tracking. Has to be a multiple of NbFramePerJob

flag_pred=1                           # To do predictive tracking. If 0 do closest neighbour tracking
flag_conf=1                           # To resolve conflict when two particles belong to the same tracjectory. Only the closest is kept

Session_INPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where tge file rays_out_ccp.hdf5 are 
Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the track_x_x.hdf5 will be  

CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Tracking3D"		# Directory where the compile file "run_submission_matlab.sh" is 
LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Tracking/LOG" 	#log directory 
OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Tracking/OUT"		#matlab output 
##############################################

export OMP_NUM_THREADS=1
for start in `seq $first $NbFramePerJobTracking $last`;
do
    stop=$(( start + NbFramePerJobTracking - 1 ));
    if [ $stop -gt $last ];
    then
        stop=$last;
    fi
    #if [ $start == 0 ];                                                                                                                                                                                                                                   
    #then                                                                                                                                                                                                                                                  
    #    continue;                                                                                                                                                                                                                                         
    #fi                                                                                                                                                                                                                                                    
    /usr/bin/qsub -cwd -V -e $LOG_path/tracking-$start-$stop-$maxdist.log -o $OUT_path/trackingOUT-$start-$stop-$maxdist.log -q piv_debian* -P PIV -N dr-$start-$stop-$maxdist  $CompileFileDir/run_submission_matlab.sh $MCRROOT "$manipname" "$NbFramePerJobMatching" "$start" "$stop" "$maxdist" "$lmin" "$flag_pred" "$npriormax" "$flag_conf" "$Session_INPUT" "$Session_OUTPUT"
done;					 
