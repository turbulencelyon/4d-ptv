#! /bin/bash                                                                                                                                                                                                                                               

####### Parameters to tune ########
ManipName="Ra1.51e10_peudense_6"		
NbFramePerJobMatching=20
Session_INPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where tge file rays.mat is 
Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the rays.dat will be saved 

CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Matching" 	#The directory where the file "submission_Centers2Rays.sh" is 
LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Centers2Rays_LOG"	#Log directory
OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Centers2Rays_OUT"	#Matlab output directory 
##########################


export OMP_NUM_THREADS=1
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/ParallelSaving_LOG.log -o $OUT_path/ParallelSaving_OUT.log -q piv_debian* -P PIV -N dr-$ManipName $CompileFileDir/run_submission_RaysSavingForParallelMatching.sh $MCRROOT "$ManipName" "$NbFramePerJobMatching" "$Session_INPUT" "$Session_OUTPUT"
	
