# GaitEventDetection 
======
The following is a MATLAB Library to detect gait events from 3-axis accelerometer signals collected during walking or running. The implemented algorithm has been published in (please remember to include an appropriate citation to acknowledge the use of library in all documents and papers that uses it):

>"Siddhartha Khandelwal; Nicholas Wickström, "**Gait Event Detection in Real-World Environment for Long-Term Applications: Incorporating Domain Knowledge into Time-Frequency Analysis**," in *IEEE Transactions on Neural Systems and Rehabilitation Engineering* , vol. 24, no. 12, pp.1363-1372, 2016"


![alt tag](http://islab.hh.se/mediawiki/images/4/4d/Figure_gaitEvents.png)


## Matlab Toolboxes required

| Name        | Version           | 
| ------------- |:-------------:| 
| MATLAB      | 9.0 | 
| Signal Processing Toolbox      | 7.2      |  
| Wavelet Toolbox | 4.16    |  
| Curve Fitting Toolbox | 3.5.3    |  

## Specifications

1. Activities: Walking and running
2. Accelerometer:
  * Placement of 3-axis Accelerometer: Anywhere around the ankle in any orientation as shown in Figure 1 [http://islab.hh.se/mediawiki/Gait_events].
  * Sensitivity of the Accelerometer: (+-) 4g or more. Please also check if the accelerometer signal has saturated during intense activity such as running.
  * Sampling Frequency: Preferred - 128 Hz [A Sampling frequency of 50Hz and above is acceptable.]
3. Input data:
  * accX - accelerometer data from X - axis
  * accY - accelerometer data from Y - axis
  * accZ - accelerometer data from Z - axis
  * Input data format: The accelerometer signals should be in units of m/s^2 and need to be in .mat format [in Matlab file format].

**IMPORTANT NOTE: The data should consist of ONLY walking and running segments of the signal. Segments corresponding to inactivity or any other activity should be removed from the signals prior to running the implementation.**

## Library details

>[HS,TO] = SK_gedAlgo(accX,accY,accZ,Fs,winSizeFactor,implement_type);

Input Arguments:

*  accX, accY, accZ [unit: m/s^2]- signals obtained from each individual axis of 3-axes accelerometer

*  Fs [unit: Hz] -   sampling frequency of the acceleration signal. Originally developed for 128 Hz

*  winSizeFactor -   size of the running window. Default value = 3. Vary this size from 2 to 6 to get better results or if code crashes.

*  implement_type -  To detect the gait event, two implementations are provided. Default value: 'fit'
  *  'fit' - 2D Gaussian distribution fitting is done to estimate gait event (runs slow) 
  *  'fast'- faster implemention to estimate the event. Warning: 'fast' implementation might not give the best estimate of the gait event compared to gaussian fitting. It is not presented in the paper but is solely provided to estimate events faster without gaussian   fitting.

Output Arguments:

*  HS - vector containing sample numbers of Heel-Strike occurences
*  TO - vector containing sample numbers of Toe-Off occurences

## MAREA Gait Database

Link to the MAREA gait database: http://islab.hh.se/mediawiki/Gait_database



