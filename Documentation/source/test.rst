Test
====

It can be usefull to test the calibration you made, because a single mistake in one plan can ruin your data processing.

The first thing to do is to visualized the point you made by using the  function ``Calib_visualisation.m``, in the Test directory, wich will overlay the point that you placed during the calibration and the raw images. If everything is fine, you should see crosses only on the black dot of the calibration plate.

An other test wich is more accurate is to generate randome particule, for wich you know the exact position, and then doing the matching process on it and see if you get the same position by matching.
The function to do this are in the test directory and the first to use is ``random_particule_generation.m``.

It will generate random particule in a large cubic volume (not the same as your mesurment volume but we will fix this later) and create a file ``Reference_point.m`` in the directory ``Calibration_test``. 

Then you have to use ``RayTracingForRandomPart.m`` to trace the rays and end up with a classical rays.mat file. Inside this file the program will consider that a plan is a frame, so you can doing the matching on it. 

So you do the matching and you plot the point obtained by this process and the reference point (a good idea is to use cross for one and circle for the second, to see the difference) and see if everything is normal.

The benefit of this :
	- you will see the gost particle 
	-you can test parameters to have an idea of the result you'll get 
	-you can see if a plan is not detect by the matching. If for example you see that for the plan 4, almost no random particles are detect by the matching process, it can be usefull to do an other calibration of this plane at least. 
	
