Experimental advices
=====================

Doing 4D-PTV experiment requires several skills that we try to detail in this section.

1. **What is the maximum area I could observed with my cameras?** On the pictures, particle has to occupied at least 3 px in diameter to be properly detected by our algorithm. Indeed, we look for spots with a gaussian shape so we need to see several luminous pixels for each particle. Besides, this criterion allows us to remove all spots caused by dust. That determines the maximum area you could observe. Let's note *d* the particle diameter. The previous criterion imposes d = 3 px on camera. So your camera resolution fixes the maximum area of observation. That is not completely exact especially for very small particles (d<100 microns) because particles reflect light making a bigger spot than their real size. But it gives a good order of magnitude.

2. Particle tracking is more efficient when frames are oversampled. Indeed, if each particle moves 10 pixels between two successives frames it will be impossible to know that it is the same particle, especially if you are tracking a lot of particle. **So how to choice a good framerate?**

    a. Estimate the maximum velocity you could expect in your setup. Let's note it *u*.
    b. Measure the spatial resolution of your cameras. We will note it :math:`a`. For instance, if you look at a 15 cm x 30 cm physical area with 1000 px x 2000 px cameras, you have a spatial resolution :math:`a = 15/1000 cm/px`.
    c. Noting :math:`f=1/T` the framerate, we compute the particle displacement :math:`D` between two successives frames: :math:`D = u dt/a`. To respect the second constraint, :math:`D<1px`. Taking :math:`D=1 px` we estimate
    
    .. math::
    
        f=1/dt = u/(aD).
        
    d. If you want to measure particle acceleration, you need to have very small noise in particle trajectories so it is even better to have D equals to a fraction of pixels. Typically taking a framerate of 3-4f is a good idea.

3. **What kind of light do I need?** As you oversampling is necessary, a powerful stable lighting is required. Its stability is important to not have issues during background soustraction. Using fluorescent particles combined with color filters could be useful because it makes possible to remove physically spots created by dust.

4. **How many particles can I track with this library?** It depends on the time you have. Less than 500 particles is very easy to track. Up to 2000 particles it is more complicated because you will need to use ``C++`` routine. Up to 10000 is more tricky because parallelisation is compulsory if you want to get trajectories within few days and not within some years!
