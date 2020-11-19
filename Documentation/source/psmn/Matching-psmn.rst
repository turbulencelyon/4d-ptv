.. _MatchingPSMN:

Matching
=========
        


How to split the matching into many jobs on the PSMN?
-----------------------------------------------------

When you try to track too many particles, run one job to do matching for the whole experiment is definetely too long. Typically, you may have some troubles beyond 3000 particles on your pictures. However the matching step can be done separatively for each frame. Indeed, we only need rays from the current frame. So it is possible to split the big job into small ones, doing matching only for some frames. That is very insteresting if you have access to a computational center which can provide you many cores simultaneously. This is what we call *parallelisation*. 

It is mot efficient (possible) for many jobs to access to a single file. But each job need for rays data to process its frames. That is why we will split rays data into smaller files. Then each job will have its own rays data file. This data spliting is realized by ``RaysSavingForParallelMatching.m`` which takes 4 arguments:

- **session**               : Path to the achitecture root,
- **ManipName**             : Name of the folder experiment,
- **camID**                 : List of cameras number,
- **NbFramePerJobMatching** : number of frames per job. Pay attention, has to be chosen as a function of processing time of one picture, in order to that each job runs for 10 min (PSMN requirements).

``RaysSavingForParallelMatching.m`` function creates a folder `Parallel/Matching/Rays/` and saves there all splitted *rays.dat* file.

.. note:: 
    With Test Data, in a matlab terminal : 
    
    .. code-block:: matlab
    
        session.path = "MyPath";
        RaysSavingForParallelMatching(session,"MyExperiment",[1,2,3],10)
        
    It will split *rays.dat* file into small files composed of 10 frames.
    
.. note:: 
    You can also launching this step in parallel at the PSMN, by compiling the function ``RaysSavingForParallelMatching.m`` and use the bash script ``submission_RaysSavingForParallelMatching.sh``. It the same way to do it that the ones describes in CenterFinding and Center2Rays. This could be usefull if you have several run to treat.
        
.. warning::
    To run your jobs on PSMN computers, it is preferable to run short jobs with a typical runtime of 10 min. The *NbFramePerJob* parameter is determined by that kind of constraints.
    
Following this method, we generate several hundreds of jobs: it it definetely not possible to run it manually. We create a ``.sh`` file which will run all jobs when it is executed. This file is created by the function ``ParallelJobsMatching.m`` which requires 11 arguments:

- **session**                : Paths to the achitecture root
- **ManipName**              : Name of the folder experiment
- **nframes**                : Total number of frames in the experiment
- **NbFramePerJobMatching**  : Number of frames per job. Pay attention, has to be chosen as a function of processing time of one picture, in order to that each job runs for 10 min (PSMN requirements).
- **CamMatch**               : Minimum number of rays to get a match
- **MaxDistance**            : Maximal authorized distance between rays to consider having a match
- **nx,ny,nz**               :  number of voxels in each direction
- **MaxMatchesPerRay**       : Maximum number of matches for one ray. 2 to consider particle overlap
- **bminx,bmaxx**            : x limits of bounding box
- **bminy,bmaxy**            : y limits of bounding box
- **bminz,bmaxz**            : z limits of bounding box
- **MinDistMatchperRay**     : Specify a volume in which you cannot have an other match if you have already found one (avoid to consider several matches for one particule)
- **Queue**                  : Running queue. By defaut it is equal to 'PIV'. It is possible to run jobs on `monointeldeb128` or `monointeldeb48` for example. Do ``qstat -g c`` to get all opened queues.

The function ``ParallelJobsMatching.m`` creates a `Parallel` folder with two subfolders `SH` and `LOG` which will contains all `sh` and `log` files for each job. The `log` file is made of all jobs output and allows you to understand what happens in case of errors. The `sh` file contains all information to run the job properly on a specific queue. This file is very specific to the PSMN.

How to choose queue?
    It is possible to see all queues and their avalaibility doing
    
    .. code-block:: bash
        
        qstat -g c

    Pay attention some queue are reserved for multi-processors jobs which is not our case. Run your jobs only on single processor queues. When you have lots of jobs, do not hesite to write to PSMN staff and ask for more cores.
    
    
The function ``ParallelJobsMatching.m`` creates also a file `<ManipName>-ParallelMatching.sh` in the folder `session.output_path/Processed_DATA/ManipName` (where ManipName is the name of the experiment). To run all jobs on the PSMN, o in a terminal

.. note:: 
    With TestData
    
    .. code-block:: bash

        cd Processed_DATA/MyExperiment/Parallel/
        sh MyExperiment-ParallelMatching.sh
    
What can I do when some jobs fail?
    It is possible to run again only these jobs doing in the SH folder:
    
    .. code-block:: bash
        
        qsub rays_n-m.sh 
        
    with n and m are the proper integers.

Matching script will save all matching files in folder `Parallel/Rays/`.
