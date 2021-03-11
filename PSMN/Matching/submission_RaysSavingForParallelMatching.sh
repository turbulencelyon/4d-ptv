#! /bin/bash                                                                                                                                                                                                                                               

####### Parameters to tune ########
ManipName="Ra3.3e10_peudense_1"		
NbFramePerJobMatching=20
Session_INPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the file rays.mat is 
Session_OUTPUT="/Xnfs/convection/Lagrangien_David_Dumont/"		#The path of the PROCESSED_DATA directory, where the rays.dat will be saved 

CompileFileDir="/Xnfs/convection/Lagrangien_David_Dumont/Toolbox/4d-ptv/Matching" 	#The directory where the file "submission_Centers2Rays.sh" is 
LOG_path="/Xnfs/convection/Lagrangien_David_Dumont/Processed_DATA/Ra3.3e10_peudense_1/LOG"	#Log directory
OUT_path="/Xnfs/convection/Lagrangien_David_Dumont/Processed_DATA/Ra3.3e10_peudense_1/OUT"	#Matlab output directory 
##########################


export OMP_NUM_THREADS=1
                                                                                                                                                                                                                                                 
/usr/bin/qsub -cwd -V -e $LOG_path/ParallelSaving_LOG.log -o $OUT_path/ParallelSaving_OUT.log -q piv_debian* -P PIV -N dr-$ManipName $CompileFileDir/run_submission_RaysSavingForParallelMatching.sh $MCRROOT "$ManipName" "$NbFramePerJobMatching" "$Session_INPUT" "$Session_OUTPUT"
	
