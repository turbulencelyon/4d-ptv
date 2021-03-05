Center finding
===============

Requirements
-------------

- Pictures of your particles from each of your cameras.

How does it work?
    This step aims to detect each particle on the pictures. The function ``CenterFinding2D.m`` will do that for you.

The particule detection is better when the background is removed. Several type of Background are calculated by ``BackgroundComputation.m``:

- BackgroundMean : it is the mean intensity value over the whole experiment in each point of pictures,
- BackgroundMax : it is the max intensity value over the whole experiment in each point of pictures,
- BackgroundMin : it is the min intensity value over the whole experiment in each point of pictures.

.. warning:: 

    In our experiment, the best detection is obtained using BackgroundMean but that can depend on your setup. We advice you to do some test with different background in order to have pictures with the clearer distinction between particles and background.
    
The function ``BackgroundComputation.m`` requires 5 arguments:

- **session**                : structure containing paths of MyPath folders,
- **ManipName**              : name of the experiment,
- **CamNum**                 : camera number,
- **StartFrame**             : number of the first frame,
- **EndFrame**               : number of the last frame,
- **Step (optional)**        : step between two frames. By defaut ajusted to compute background over 1000 frames uniformly reparted in the experiment,
- **format (optional)**      : pictures name format : equals to '%05d.tif' by defaut. The beginning of picture names has to be '%ManipName_cam%CamNum_%format'.

It is not necessary to use all frames to compute background. Indeed, we oversample so successive frames are highly correlated. Averaging over decorrelated frames is sufficient to get a good background in a reasonable delay.

The function ``BackgroundComputation`` saves the three backgrounds in the file *session.output_path/Processed_DATA/%ManipName/Background_camcam%NumCam.mat*.

.. Warning:: 

    **session.input_path** and **session.output_path** arguments define . In a basic installation, both ot them are equal to **MyPath/** but it could be different if your *DATA/* and *Processed_DATA/* directories are not located in the same *MyPath/* folder.

.. note::

    To compute background of test data, try in a matlab terminal

    .. code-block:: matlab
        
        session.input_path = "My4DPTVInstallationPath/Documentation/TestData/";  % My4DPTVInstallationPath has to be adapted !!!
        session.output_path = "My4DPTVInstallationPath/Documentation/TestData/";
        BackgroundComputation(session,"MyExperiment",1,1,100,10) # for camera 1
        BackgroundComputation(session,"MyExperiment",2,1,100,10) # for camera 2
        BackgroundComputation(session,"MyExperiment",3,1,100,10) # for camera 3

Once the background is calculated, you can launch ``CenterFinding2D.m``. This function will treat all pictures for a camera. You have to launch it for each camera. ``CenterFinding2D.m`` takes at least 6 arguments:

- **session**                   : structure containing paths of MyPath folders,
- **ManipName**                 : name of the experiment,
- **NumCam**                    : camera number,
- **firstFrame**                : the first frame to treat,
- **nframes**                   : the last frames to treat,
- **th**                        : threshold value to detect points,
- **sz**                        : typical point diameter (in pixels),
- **Test (optional)**           : true/false. If true enters in test mode,
- **BackgroundType (optional)** : determine which background is substracted to pictures. By defaut is equal to BackgroundMean,
- **format (optional)**         : pictures name format : equals to '%05d.tif' by defaut. The beginning of picture names has to be '%ManipName_cam%CamNum_%format'. 

The threshold and point diameter values depend on the camera treated. To set the best values, a test mode can be activated thanks to *Test* argument. In test mode, only the first frame is analysed and several plots are shown to help you to determine the best values. We advice you to set these values in order to get more or less the same number of detected particles for each camera.

``CenterFinding2D.m`` creates a file *session.output_path/Processed_DATA/%ManipName/centers_cam%NumCam.mat* which contains all detected points.



.. note::

    To use testData, test the following lines in a matlab terminal:

    .. code-block:: matlab
    
        session.input_path = "My4DPTVInstallationPath/Documentation/TestData/";  % My4DPTVInstallationPath has to be adapted !!!
        session.output_path = "My4DPTVInstallationPath/Documentation/TestData/";
        CC1 = CenterFinding2D(session,"MyExperiment",1,1,100,6000,3) # for camera 1
        CC2 = CenterFinding2D(session,"MyExperiment",2,1,100,6000,3) # for camera 2
        CC3 = CenterFinding2D(session,"MyExperiment",3,1,100,6000,3) # for camera 3


.. warning::

    This step can be experiment-dependant. For your particular experiment, what we provide for this step could be not efficient. However, the following steps are more general.

How to run a compiled version of the ``CenterFinding2D.m``?
-------------------------------------------------------------

It is possible to compile ``CenterFinding2D.m`` function to run it outside a MATLAB instance directly in a terminal. This can be useful to run it on cluster, for instance at the PSMN. The matlab function to use for that is ``submission_center_finding.m``. 

    1. If you don't have the compiled files yet (an executable ``submission_center_finding`` and a bash script ``run_submission_center_finding.sh``), compile the script ``submission_center_finding.m`` doing in a matlab terminal:


        .. code-block:: matlab
            
            mcc -m submission_center_finding.m
            
        An executable file ``submission_center_finding`` and a bash file ``run_submission_center_finding.sh`` will appear in the same folder.

    2. Modifie the line 30 of the file "run_submission_center_finding" to add the path of the executable file like this:
    
        .. code-block:: bash

              eval "/MyPath/Submision_center_finding" $args
              
    
    3. To run it in your machine:

        .. code-block:: bash

            sh run_submission_center_finding.sh $MCRROOT "ManipName" "CamNum" "FirstFrame" "Nframes" "Th" "Size" "Session_INPUT" "Session_OUTPUT"
            
.. warning:: 
    
    Even if some parameters are numbers (integers or floats), you need to tipe them as string by using the quote ".
    
    
If you want to launch the function at the PSMN and use parallelisation, use the file ``submission_CenterFinding.sh``:
    
    1. Change the parameters at the begining of the script to use your own parameters 
    
        .. code-block:: bash            
            
            ManipName="MyExperiment"    
            CamNum=3                                #The camera on wich you want to find the center 
            FirstFrame=300                          #The first frame (usefull if you don't start at one)
            Nframes=36000                           #The final frame to treat 
            Th=6500                                 #Threshold to detec a part (it have to be tune with the function CenterFinding.m and with test=true 
            Size=5                                  #The size of a part (in pixel)
            Session_INPUT="/MyWorkspace/"       #The path of the DATA directory, where all the image are 
            Session_OUTPUT="/MyWorkspace/"      #The path of the PROCESSED_DATA directory, where the centercamk.mat will be saved

            CompileFileDir="/home/eberna07/Stage_EB_2020/4d-ptv/CenterFinding"          #The directory where the file "runSubmision_center_finding.sh" is 
            LOG_path="/MyWorkspace/MyExperiment/CenterFinding_LOG"      #log directory (warning: the directory has to be created before launch the code)
            OUT_path="/MyWorkspace/MyExperiment/CenterFinding_OUT"      #matlab output (warning: the directory has to be created before launch the code)

        
    2. Launch this function in a terminal doing:
  
        .. code-block:: bash
            
            sh submission_CenterFinding.sh  
            
        This will launch a job at the PSMN, on the queue PIV, you can check if everything is ok by looking at the file ``center_camCamNum.log`` in the LOG directory.
        
.. warning:: 

    To use PSMN installations see :ref:`CenterFindingPSMN`

