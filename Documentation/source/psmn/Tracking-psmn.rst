.. _TrackingPSMN:

Tracking
==========


How to split Tracking into many jobs on the PSMN?
-------------------------------------------------

As for the matching, as soon as you track several thousands of particles, it becomes impossible to get trajectories in a raisonnable time (reasonable meaning a few days). However, the tracking is harder to parallelize because it requires matched points for all frames. So instead of running one job for tracking, we run several jobs tracking particles only over several hundreds of frames. The precise number of frames treated by each job is, as previously, determined by the required time to do it.
Doing that, we have particules trajectory for each sets of frames. To get the complete trajectories over the whole experiment, there need to be reconnected and that is done in the following step: the *stitching*.

The function ``track3d_psmn.m`` computes trajectory like previously but it loads several matched files simultaneously. Indeed, the tracking is faster than matching so we can treat more frames per job. Besides, as it is not possible to estimate next particle position when we have its position for less than *mpriormax* frames. In that case, we do simple closest neighbour tracking so  it makes no sense to track a small number of frame per job with this code because otherwise we do not do predictive tracking but we do closest neighbour tracking. That's why we load several matched files. The function ``track3d_psmn.m`` takes 10 arguments:

- **session**               : Path to the achitecture root
- **ManipName**             : Name of the experiment
- **NbFramePerJobTracking** : Number of frame per job for parallelized matching
- **minframe**              : number of the first frame to treat
- **maxframe**              : number of the last frame to treat. Pay attention, min and max frame have to be multiple of NbFramePerJob
- **maxdist**               : maximum travelled distance between two successive frames
- **lmin**                  : minimum length of a trajectory (number of frames)
- **flag_pred**             : 1 for predictive tracking, 0 otherwise
- **npriormax**             : maximum number of prior frames used for predictive tracking
- **flag_conf**             : 1 for conflict solving, 0 otherwise

.. figure:: Figures/InOutputtrack3d_psmn.png
    :width: 110%
    
    Input and output files of ``track3d_psmn.m`` function.

As previously, trajectories are saved in a file */Parallel/Tracking/Tracks/tracks_{minframe}-{maxframe}.h5*. This file can also be read with `readmatches.m` function.

It is better to compile ``track3d_psmn.m`` function.
	1. Again, if you don't have the compiled file yet, compile the function ``submission_Tracking.m``
	
		.. code-block:: matlab
			
			mcc -m submission_Tracking.m
			
	2. Modify the line 30 of the file ``run_submission_Tracking.sh`` to add the path of the executable file like this:
	
        .. code-block:: bash

              eval "/MyPath/submission_Tracking" $args
              
    3. To run it in your machine:

        .. code-block:: bash

            sh run_submission_Tracking.sh "MyExperiment" "NbFramePerJobMatching" "FirstFrame" "LastFrame" "maxdist" "lmin" "flag_pred" "npriormax" "flag_conf" "Session_INPUT" "Session_OUTPUT"

.. warning:: 
	
	Even if some parameters are numbers (integers or floats), you need to tipe them as string by using the quote ".
	
	
To run all jobs simultaneously use ``submission_Tracking.sh`` file after completing its header:

		.. code-block:: bash
		
			NbFramePerJobMatching=20             # Number of frame per job for parallel matching
			maxdist=0.4                         # maximum distance between rays to consider a match
			lmin=5                                # minimum trajectory length
			npriormax=5                           # number of points used to predict next particle position
			manipname="Ra1.51e10_peudense_6"
			first=401                               # First frame of the experiment                                                                                           
			last=36000                           # Last frame of the experiment
			NbFramePerJobTracking=5000             # Number of frame per job for tracking. Has to be a multiple of NbFramePerJob

			flag_pred=1                           # To do predictive tracking. If 0 do closest neighbour tracking
			flag_conf=1                           # To resolve conflict when two particles belong to the same tracjectory. Only the closest is kept

			Session_INPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where tge file rays_out_ccp.hdf5 are 
			Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the track_x_x.hdf5 will be  

			CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Tracking3D"		# Directory where the compile file "run_submission_matlab.sh" is 
			LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Tracking/LOG" 	#log directory 
			OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Tracking/OUT"		#matlab output 



	Several parameters are very important:

			- **minframe** and **maxframe** the first and last are number of the first and last frames of the experiment,
			- **NbFramePerJobMatching** is the number of frame per job for parallel matching,
			- **NbFramePerJobTracking** is the number of frame per job for parallel tracking: it has to be a multiple of **NbFramePerJobMatching** because it will open several matching output files until achieves **NbFramePerJobTracking**. This number has to be selected as a function of computational time. Typically it is equal to several thousands.

	Once the submission files completed, you can launch it by opening a terminal in the the tracking directory and tipeing the command     
    
        .. code-block::
        
            sh submission_Tracking.sh 

.. note::
    You can see if your job are running by doing ``qstat``. If their state are ``qw`` it mean that all the CPU of the queue are running and your job is in waiting state. Then, if everything is ok, you will see the state ``r``, meaning that the job is running. If you see ``eqw``, it means that there is a problem but you can info on this problem by tiping the command ``qstat -explain E -j`` and the number of the job. In general, it's because the log and out directory you have defined are not created.
    The exact path depends on where you are precisely in the folder. We precise that is is not necessary to parallelize tracking for test data as data are very small, it is presented only to understand processing.

