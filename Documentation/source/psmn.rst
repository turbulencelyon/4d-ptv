****************************************************************
How to use PSMN installations to accelerate the data processing?
****************************************************************
At the ENS de Lyon, we have access to PSMN installations which gives us access to several thousands of cores. As soon as it is possible, we can run many short jobs simultaneously to accelerate data processing. Our codes were made for that and the following documentation will show you how does it work. All configuration files are specific to PSMN but our codes can be used everywhere with configuration corresponding to your specific computational system.

.. warning::
    All functions specific to the PSMN are in the folder ``PSMN`` of the module. In this folder, you will find a
    folder for each step of the data processing (CentersFinding, Centers)
    Rays...) and in these folders, the specific functions for the step.

.. include:: psmn/CentersFinding-psmn.rst

.. include:: psmn/Centers2Rays-psmn.rst

.. include:: psmn/Matching-psmn.rst

.. include:: psmn/Tracking-psmn.rst

.. include:: psmn/Stitching-psmn.rst
