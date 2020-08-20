#! /bin/bash 

                                                                                                                                                                                                                                              
############  Parameters to tune  ######################
ManipName="Ra1.51e10_peudense_2"	
CamNum=3								#The camera on wich you want to find the center 
FirstFrame=300							#The first frame (usefull if you don't start at one)
Nframes=36000							#The final frame to treat 
Th=6500						#Threshold to detec a part (it have to be tune with the function CenterFinding.m and with test=true 
Size=5									#The size of a part (in pixel)
Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the DATA directory, where all the image are 
Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the centercamk.mat will be saved

CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/CenterFinding"			#The directory where the file "runSubmision_center_finding.sh" is 
LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_2/Center_Finding_LOG"	#log directory
OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_2/Center_Finding_OUT"	#matlab output 

##############################################
export OMP_NUM_THREADS=1
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/center_cam$CamNum.log -o $OUT_path/center_cam$CamNum.log -q piv_debian* -P PIV -N dr-cam$CamNum $CompileFileDir/run_Submission_center_finding.sh $MCRROOT "$ManipName" "$CamNum" "$FirstFrame" "$Nframes" "$Th" "$Size" "$Session_INPUT" "$Session_OUTPUT"
	
