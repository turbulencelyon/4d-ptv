#! /bin/bash                                                                                                                                                                                                                                               

############  Parameter to tune  ######################
ManipName="Ra1.51e10_peudense_TOTAL"
SeuilInit=1
DeltaSeuil=0.0001
DeltaInc=0.0001
Epsilon=0.00000464
Ech=150
NbFrame=200
MinConv=100

Session_INPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where tge file rays_out_ccp.hdf5 are 
Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the track_x_x.hdf5 will be  

CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Post-Processing/Pair_dispersion_EBernard/"		# Directory where the compile file "run_submission_matlab.sh" is 
LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_TOTAL/Post_Processed_Data/PairDipsersion/MinSize_450/PairDispFinal/LOG" 	#log directory 
OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_TOTAL/Post_Processed_Data/PairDipsersion/MinSize_450/PairDispFinal/OUT"		#matlab output 
##############################################

export OMP_NUM_THREADS=1
                                                                                                                                                                                                                                                  
/usr/bin/qsub -cwd -V -e $LOG_path/PariDisp-$SeuilInit.log -o $OUT_path/PariDisp-$SeuilInit.log -q piv_debian* -P PIV -N PairDisp-$SeuilInit  $CompileFileDir/run_submission_PairDisp.sh $MCRROOT "$Session_INPUT" "$Session_OUTPUT" "$ManipName" "$SeuilInit" "$DeltaSeuil" "$DeltaInc" "$Epsilon" "$Ech" "$NbFrame" "$MinConv" 
