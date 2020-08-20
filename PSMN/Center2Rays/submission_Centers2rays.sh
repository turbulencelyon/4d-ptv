#! /bin/bash                                                                                                                                                                                                                                               

####### Parameters to tune ########
kcam=4		#The camera on wich you want to trace the rays
CalibPath="/Xnfs/convection/Stage_EB_2020/New_Calibration/calib.mat"	#The path of the calibration file
ManipName="Ra1.51e10_peudense_6"		
FirstFrame=400
Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the "center_camX.mat" are 
Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the "rays_camX.mat" will be saved

CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Center2Rays" 	#The directory where the file "submission_Centers2Rays.sh" is 
LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Centers2Rays_LOG"	#Log directory
OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Centers2Rays_OUT"	#Matlab output directory 
##########################


export OMP_NUM_THREADS=1
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/center2rays_$kcam.log -o $OUT_path/center2rays_OUT$kcam.log -q piv_debian* -P PIV -N dr-center2rays-$kcam $CompileFileDir/run_submission_Centers2Rays.sh $MCRROOT "$kcam" "$CalibPath" "$ManipName" "$FirstFrame" "$Session_INPUT" "$Session_OUTPUT"
	
