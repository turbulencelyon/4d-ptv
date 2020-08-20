Post-Processing
===============

Once you have created the track, we can use them to obtain physical result. 

Pair dispersion 
---------------

One interesting thing to look on while doing 4d-PTV is relative pair dispersion. It's easy to do since we have the position of the particules at any time. 
	
	
.. figure:: Figures/PairDisp.png
    :width: 60%
	
	
	
To compute pair dispersion, you'll need three step. 

First, you ned to compare the initial separation of all the couple of particules. To do this, use the function ``Dinitial.m``, wich take three argument
		
- **session** 		: the path of your Processed_DATA directory. 
- **ManipName** 	: the name of the experimentation 
- **size** 			: the minimum under wich you will not consider a track. 
- **track** 		: the structure that contains on wich you want to compute the initial separation
		
This step will provide you 3 files: ``Dinitial.mat`` wich contains the initial separation between all the track, ``S.mat`` and ``index.mat`` that contains the position in the S structure on wich we are interested.
	
Then, you can use the function ``PariDisp.m`` wich will compute the pair dispersion for a given inital separation. The arguments are:
		
- **session** 	: the path of your Processed_DATA directory. 
- **ManipName** 	: the name of the experimentation 
- **SeuilInit** 	: the initial separation you want to compute
- **DeltaSeuil** : the minimum size of the bin. This parameters can varie in function of the number of track you want for convergence 
- **DeltaInc** 	: the increment you want on the bin size if the initial one is too small to get convergence criteria (typically DeltaSeuil/2 is a good value).
- **Espilon** 	: the rate of kinetic energy dissipation in your experiment 
- **Ech** 		: the acquisition rate (in Hz)
- **MinConv** 	: the minimum of couple of track per bin. This will have an impact on the size of the bin, if you aksed for too much particles, your bin size could become too large
- **NbFrame** 	: the value of the ``size`` argument you used in ``Dinitial.mat`` 
		
At the end of this step, you will obtain a structure containing the statistic of pair dispersion, the bin size, the time and the value of t* and the number of couple of tracks as a function of time. This last data is very important because the number of couple will decrease with time (because you will have less and less of long track) and this could strongly destroy the convergence of your curve. You will see it when ploting the pair dispersion, at first it's very smooth curve but as you loose tracks at every time step, the curve will become noisy and small jump are going to appear, even that after a certain time the data mean nothing. So pay attention to this when you interpret your curve.
You can also use a compiled version of the function using:
	
	.. code-block:: matlab 
	
        mcc -m submission_PairDisp.m	

Then, it possible to launch job at the PSMN, as for the previous step, by using ``submission_PairDisp.sh`` and complete the header of the function. 
	
Once you've done the previous processing step, using can use the function ``PlotDispersion.m`` to plot the pair dispersion. It takes as arguments:
		
- **session**		: the path of your Processed_DATA directory. 
- **ManipName** 		: the name of the experimentation 
- **Init** 			: the first initial separation to plot
- **Last** 			: the last initial separation to plot 
- **Pas** 			: the step between 2 initial separation 
- **Norm** (optional): it's a boolean (``true`` of ``false``) to set to ``true`` (it's ``false`` if there is no indication) if you want to compensate the curbe by t^n
- **P**  			: use it if you set ``Norm`` to ``true``, P will be the exponant of t.
