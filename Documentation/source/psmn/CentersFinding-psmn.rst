.. _CenterFindingPSMN:

Center finding
==============


How to run a compiled version of the ``CenterFinding2D.m``?
-----------------------------------------------------------

It is possible to compile ``CenterFinding2D.m`` function to run it outside a MATLAB instance directly in a terminal. This can be useful to run it on cluster, for instance at the PSMN. The matlab function to use for that is ``submission_center_finding.m``. 

    1. If you don't have the compiled files yet (an executable ``submission_center_finding`` and a bash script ``run_submission_center_finding.sh``), compile the script ``submission_center_finding.m`` doing in a matlab terminal:


        .. code-block:: matlab
            
            mcc -m submission_center_finding.m
            
        An executable file ``submission_center_finding`` and a bash file ``run_submission_center_finding.sh`` will appear in the same folder.

	2. Modifie the line 30 of the file "run_submission_center_finding" to add the path of the executable file like this:
	
        .. code-block:: bash

              eval "/MyPath/Submission_center_finding" $args
              
	
    3. To run it in your machine:

        .. code-block:: bash

            sh run_submission_center_finding.sh $MCRROOT "ManipName" "CamNum" "FirstFrame" "Nframes" "Th" "Size" "Session_INPUT" "Session_OUTPUT"
            
.. warning:: 
	
	Even if some parameters are numbers (integers or floats), you need to tipe them as string by using the quote ".
	
How to split CentersFinding into many jobs on the PSMN?
-------------------------------------------------------

If you want to run the function at the PSMN and use parallelisation, use the file ``submission_CenterFinding.sh``:
	
	1. Change the parameters at the begining of the script to use your own parameters 
	
        .. code-block:: bash            
            
			ManipName="MyExperiment"	
			CamNum=3								#The camera on which you want to find the center 
			FirstFrame=300							#The first frame (useful if you don't start at one)
			Nframes=36000							#The final frame to treat 
			Th=6500									#Threshold to detect a part (it has to be tuned with the function ``CenterFinding.m`` and with test=true 
			Size=5									#The size of a part (in pixel)
			Session_INPUT="/MyWorkspace/"		#The path of the DATA directory, where all the images are 
			Session_OUTPUT="/MyWorkspace/"		#The path of the PROCESSED_DATA directory, where the centercamk.mat will be saved

			CompileFileDir="MyPath/4d-ptv/CenterFinding"			    #The directory where the file "runSubmision_center_finding.sh" is 
			LOG_path="/MyWorkspace/MyExperiment/CenterFinding_LOG"		#log directory (warning: the directory has to be created before launch the code)
			OUT_path="/MyWorkspace/MyExperiment/CenterFinding_OUT"		#matlab output (warning: the directory has to be created before launch the code)

        
	2. Run this function in a terminal doing:
  
        .. code-block:: bash
            
            sh submission_CenterFinding.sh  
            
        This will launch a job at the PSMN, on the queue PIV, you can check if everything is ok by looking at the file ``center_camCamNum.log`` in the LOG directory.

