#! /bin/bash                                                                                                                                                                                                                                               

####### Parameters to tune ########
Maxkcam=4		#The camera on which you want to draw the rays
CalibPath="/Xnfs/convection/Stage_EB_2020/New_Calibration/calib.mat"	#The path of the calibration file
ManipName="Ra3.3e10_peudense_2"		
FirstFrame=271
Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the "center_camX.mat" are 
Session_OUTPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the "rays_camX.mat" will be saved

CompileFileDir="/Xnfs/convection/Lagrangien_David_Dumont/Toolbox/4d-ptv/Center2Rays" 	#The directory where the file "submission_Centers2Rays.sh" is 
LOG_path="$Session_OUTPUT/Processed_DATA/$ManipName/LOG"	#Log directory
OUT_path="$Session_OUTPUT/Processed_DATA/$ManipName/OUT"	#Matlab output directory 
##########################


export OMP_NUM_THREADS=1
for kcam in `seq 1 $Maxkcam`;
do                                                                                                                                                                                                                        
/usr/bin/qsub -cwd -V -e $LOG_path/center2rays_$kcam.log -o $OUT_path/center2rays_OUT$kcam.log -q CLG5218deb192*,CLG6226Rdeb192* -N C2R-$kcam-$ManipName $CompileFileDir/run_submission_Centers2Rays.sh $MCRROOT "$kcam" "$CalibPath" "$ManipName" "$FirstFrame" "$Session_INPUT" "$Session_OUTPUT"
done	
