#! /bin/bash                                                                                                                                                                                                                                               

############  Parameter to tune  ######################
NbFramePerJobMatching=150             # Number of frame per job for parallel matching
maxdist=0.7                         # maximum distance between rays to consider a match
lmin=5                                # minimum trajectory length
npriormax=5                           # number of points used to predict next particle position
manipname="Ra4.6e9_peudense2_6"    #"Ra1.51e10_peudense_7"
first=1                              # First frame of the experiment                                                                                           
last=36000                          # Last frame of the experiment
NbFramePerJobTracking=6000             # Number of frame per job for tracking. Has to be a multiple of NbFramePerJob

flag_pred=1                           # To do predictive tracking. If 0 do closest neighbour tracking
flag_conf=1                           # To resolve conflict when two particles belong to the same tracjectory. Only the closest is kept

Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the files rays_out_ccp.hdf5 are 
Session_OUTPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the track_x_x.hdf5 will be  

CompileFileDir="/Xnfs/convection/Lagrangien_David_Dumont/Toolbox/4d-ptv/PSMN/Tracking3D"		# Directory where the compiled file "run_submission_matlab.sh" is 
LOG_path="$Session_OUTPUT/Processed_DATA/$manipname/Parallel/Tracking/LOG" 	#log directory 
OUT_path="$Session_OUTPUT/Processed_DATA/$manipname/Parallel/Tracking/OUT"		#matlab output 
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
    /usr/bin/qsub -cwd -V -e $LOG_path/tracking-$start-$stop-$maxdist.log -o $OUT_path/trackingOUT-$start-$stop-$maxdist.log -q CLG5218deb192*,CLG6226Rdeb192* -N track-$start-$stop-$maxdist  $CompileFileDir/run_submission_Tracking.sh "/applis/PSMN/generic/Matlab/R2019b" "$manipname" "$NbFramePerJobMatching" "$start" "$stop" "$maxdist" "$lmin" "$flag_pred" "$npriormax" "$flag_conf" "$Session_INPUT" "$Session_OUTPUT"
done;					 
