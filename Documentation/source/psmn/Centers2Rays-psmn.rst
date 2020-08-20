.. _Centers2RaysPSMN:

Center to rays
================

How to parallelize this step?
-----------------------------------

	It's possible to split the work for each camera instead of doing all the camera in a row, for this you can use the function ``Centers2RaysParallel`` which will give a number x file center_camx.mat, with x the number of camera. 

	Once you have those file for all the camera, you can use the function ``Ray_recombinaison.m`` which will gathered all the rays in one file ``rays.mat``
	
``Center2RaysParallel`` takes 6 arguments:  

- **session**            : structure containing paths of MyPath folders,
- **ManipName**          : name of the experiment,
- **Calib**              : calib.mat file,
- **kcam**               : the number of the cam you want to treat,
- **FirstFrame**		 : the first frame to treat 
- **Ttype** (optional)   : type of the transformation to use. 'T1' for linear transformation (defaut). 'T3' for cubic transformation.


``Rays_recombinaison.m`` takes 3 arguments:

- **session**            : structure containing paths of MyPath folders,
- **ManipName**          : name of the experiment,
- **camID**              : list of camera numbers. ex: [1,2,3] if you have 3 cameras numbered 1,2,3 respectively,


How to run a compiled version of  ``Center2Rays``?
--------------------------------------------------

It can be useful to run this step at the PSMN and for this you can compile the Matlab function and use the bash submission function. The matlab function to use is ``submission_Centers2Rays.m``.

	1. If you don't have the compiled files yet (an executable ``submission_Centers2Rays`` and a bash script ``run_submission_Centers2Rays``), compile the script ``submission_Centers2Rays.m`` doing in a matlab terminal:


		.. code-block:: matlab
			
			mcc -m submission_Centers2Rays.m -a /applis/PSMN/generic/Matlab/R2017b/toolbox/images/images
			
	2. Modify the line 30 of the file ``run_submission_center_finding.sh`` to add the path of the executable file like this:
	
        .. code-block:: bash

              eval "/MyPath/submission_Centers2Rays" $args
              
    3. To run it in your machine:

        .. code-block:: bash

            sh run_submission_center_finding.sh $MCRROOT "$kcam" "$CalibPath" "$ManipName" "$FirstFrame" "$Session_INPUT" "$Session_OUTPUT"

How to split this step into many jobs on the PSMN?
--------------------------------------------------

If you want to run the function at the PSMN and use parallelisation, once you have executed the previous paragraph, use the file ``submission_Centers2Rays.sh``:
	
	1. Change the parameters at the begining of the script to use your own parameters
	 
        .. code-block:: bash
        
			kcam=4		#The camera on which you want to trace the rays
			CalibPath="/MyPath/calib.mat"	#The path of the calibration file
			ManipName="MyExperiment"		
			FirstFrame=400
			Session_INPUT="/MyWorkspace/"		#The path of the PROCESSED_DATA directory, where the "center_camX.mat" are 
			Session_OUTPUT="/MyWorkspace/"		#The path of the PROCESSED_DATA directory, where the "rays_camX.mat" will be saved

			CompileFileDir="/MyWorkspace/4d-ptv/Center2Rays" 	#The directory where the file "submission_Centers2Rays.sh" is 
			LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Centers2Rays_LOG"	#Log directory
			OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Centers2Rays_OUT"	#Matlab output directory 


	2. Run this function for each camera in a terminal doing:
  
        .. code-block:: bash
            
            sh submission_Centers2Rays.sh
            
	3. Once you have all your file rays_camX.mat, launch the function ``Rays_recombinaison.m`` in a Matlab terminal doing:

		.. code-block:: matlab
			
			Rays_recombinaison(session,'MyExperiment',[1 2 3 4])
