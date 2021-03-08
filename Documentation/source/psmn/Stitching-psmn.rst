.. _StitchingPSMN:

Stitching
==========

How to split Stitching processing into many jobs on the PSMN?
--------------------------------------------------------------

When the tracking has been made by several jobs, it is necessary to reconnect trajectories in between successive files and then we can do *classic* stitching. There are one or two steps depending on the total number of frame of your experiment. The stitching is quitte time consuming so if you have many frames (more than 20000) for your experiment, it is clever to do parallel stitching on a reduced number of frames. Doing that, you will get several packets of reconnected trajectories. The second step allows you to reconnect these packets into only one. If you have few frames in your experiment (less than 10000) it is worth to run only one stitching job working on all frames.

If you have 50000 frames per experiment and that the tracking was done on
packets of 1000 frames, you can run 25 jobs doing stitching on packets of 2000
frames using ```Stitching_psmnA.m` function and then reconnect all the packet
into one using ```Stitching_psmnB.m```. 

The first step is realised by ``stitchTracksSides.m``. But the user has just to use ``Stitching_psmnA.m`` function which does everything for him. This function takes 11 arguments:

- **session**                : session.path contains MyPath,
- **ManipName**              : Name of the experiment,
- **minframe**               : First frame to process,
- **maxframe**               : Last frame to process,
- **NbFramePerJobTracking**  : Number of frame per tracking job,
- **FileName**               : Name of the tracks file,
- **dfmax**                  : maximum number of tolerated missing frames to reconnect to trajectories,
- **dxmax**                  : maximum tolerated distance (in norm) between projected point after the first trajectory and the real beginning position of the stitched one,
- **dvmax**                  : maximum tolerated relative velocity difference between the end of the first trajectory and the beginning of the stitched one,
- **lmin**                   : minimum length for a trajectory to be stitched. 

```Stitching_psmnA``` reconnectes trajectories in a small packet of a few thousands of frames. ```Stitching_psmnB``` reconnectes trajectories between the small packets: to get the trajectories from 2 successive packets, it looks for the last *dfmax* frame of the first packet and the first *dfmax* frames of the second packet and reconnect trajectories only within these frames. Indeed, trajectories were already reconnected anywhere else in the packets by ```Stitching_psmnA```.
`

This final step can be splitted into several jobs at the PSMN, in order to save time. To do it you have to:

    1. If you don't have the compiled files yet you can compile the function ``submission_Matlab_Stitching.m``

        .. code-block:: matlab
            
            mcc -m submission_Matlab_Stitching
        
    2. Modify the line 30 of the file "run_Submission_Matlab_Stitching.sh" to add the path of the executable file like this:

        .. code-block:: bash

            eval "/MyPath/Submission_Matlab_Stitching" $args

    3. To run it in your machine:

        .. code-block:: bash

            sh run_Submission_Matlab_Stitching.sh $MCRROOT "ManipName" "FirstFrame" "LastFrame" "dfmax" "dxmax" "dvmax" "lmin" "Session_INPUT" "Session_OUTPUT"
            
            
.. warning:: 
	
	Even if some parameters are numbers (integers or floats), you need to type them as string by using the quote ".
	
	
As for the Tracking, CenterFinding and Centers2Rays, you can launch several job at the PSMN by using ``submission_Stitching.sh``:
		
		1. Complete the header of the function to tune your parameters 
		
		    .. code-block:: bash 
		
			    ManipName="Ra1.51e10_peudense_6"
			    dfmax=60              # Number of frame per job for parallel matching
			    dxmax=0.5                           # maximum distance between rays to consider a match
			    dvmax=0.35                               # minimum trajectory length
			    lmin=5                           # number of points used to predict next particle position
			    first=35401                               # First frame of the experiment                                                                                           
			    last=36000                              # Last frame of the experiment
			    NbFramePerJobTracking=600              # Number of frame per job for tracking. Has to be a multiple of NbFramePerJob

			    Session_INPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where tge file track_x_x.h5 are 
			    Session_OUTPUT="/Xnfs/convection/Stage_EB_2020/"		#The path of the PROCESSED_DATA directory, where the StitchA_x_x.h5 will be  

			    CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/Stitching" #Directory where the file "run_Submission_Matlab_Stitching.sh" is
			    LOG_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Stitching/LOG"	#Log directory
			    OUT_path="/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.51e10_peudense_6/Parallel/Stitching/OUT"	#Matlab output directory

		2. Open a terminal in the directiry Stitching and use the command: 
		
		    .. code-block:: bash 
		
			    sh submission_Stitching.sh 
			
