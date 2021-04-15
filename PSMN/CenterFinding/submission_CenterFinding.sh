#! /bin/bash 
                                                                                                                                                                                                                                      
############  Parameters to tune  ######################
ManipName="Ra3.3e10_peudense_2"	
MaxCamNum=4								#The camera on wich you want to find the center 
FirstFrame=271							#The first frame (usefull if you don't start at one)
Nframes=36000							#The final frame to treat 
Th=5000						#Threshold to detec a part (it have to be tune with the function CenterFinding.m and with test=true 
Size=3									#The size of a part (in pixel)
Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont"		#The path of the DATA directory, where all the image are 
Session_OUTPUT="/Xnfs/convection/Lagrangien_David_Dumont"		#The path of the PROCESSED_DATA directory, where the centercamk.mat will be saved

CompileFileDir="/Xnfs/convection/Lagrangien_David_Dumont/Toolbox/4d-ptv/CenterFinding"			#The directory where the file "runSubmision_center_finding.sh" is 
LOG_path="$Session_OUTPUT/Processed_DATA/$ManipName/LOG"	#log directory
OUT_path="$Session_OUTPUT/Processed_DATA/$ManipName/OUT"	#output directory

##############################################
export OMP_NUM_THREADS=1
for CamNum in `seq 1 $MaxCamNum`;
do                                                                                                                                                                                                                              
/usr/bin/qsub -cwd -V -e $LOG_path/center_cam$CamNum.log -o $OUT_path/center_cam$CamNum.log -q piv_debian* -P PIV -N CF-cam$CamNum-$ManipName $CompileFileDir/run_submission_CenterFinding2D.sh $MCRROOT "$ManipName" "$CamNum" "$FirstFrame" "$Nframes" "$Th" "$Size" "$Session_INPUT" "$Session_OUTPUT"
done
	
