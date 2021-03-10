#! /bin/bash                                                                                                                                                                                                                                               

####### Parameters to tune ########
kcam=2		#The camera on which you want to trace the rays
CalibPath="/Xnfs/convection/Lagrangien_David_Dumont/DATA/CalibrationTilte2/calib.mat" #"/Xnfs/convection/Stage_EB_2020/New_Calibration/calib.mat"	#The path of the calibration file
ManipName="Ra3.3e10_peudense_1"		
FirstFrame=1
Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the "center_camX.mat" are 
Session_OUTPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the "rays_camX.mat" will be saved

CompileFileDir="/Xnfs/convection/Lagrangien_David_Dumont/Toolbox/4d-ptv/Center2Rays" 	#The directory where the file "submission_Centers2Rays.sh" is 
LOG_path="/Xnfs/convection/Lagrangien_David_Dumont/Processed_DATA/Ra3.3e10_peudense_1/LOG"	#Log directory
OUT_path="/Xnfs/convection/Lagrangien_David_Dumont/Processed_DATA/Ra3.3e10_peudense_1/OUT"	#Matlab output directory 
##########################


export OMP_NUM_THREADS=1
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/center2rays_$kcam.log -o $OUT_path/center2rays_OUT$kcam.log -q piv_debian* -P PIV -N dr-center2rays-$kcam $CompileFileDir/run_submission_Centers2Rays.sh $MCRROOT "$kcam" "$CalibPath" "$ManipName" "$FirstFrame" "$Session_INPUT" "$Session_OUTPUT"
	
