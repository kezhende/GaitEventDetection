%% 
% Copyright (c) <2016>, <Siddhartha Khandelwal>

% Publication to cite: 
% Siddhartha Khandelwal, Nicholas Wickström (2016). 
% Gait Event Detection in Real-World Environment for Long-Term Applications : 
% Incorporating Domain Knowledge into Time-Frequency Analysis. 
% IEEE transactions on neural systems and rehabilitation engineering.
%% Example main file

load accX 
load accY 
load accZ

% Input arguments

% accX, accY, accZ [unit: m/s^2]- signals obtained from each individual axis
%                                 of 3-axes accelerometer
% Fs [unit: Hz] -   sampling frequency of the acceleration signal.
%                   Originally developed for 128 Hz
% winSizeFactor -   size of the running window. Default value = 3. 
%                   Vary this size from 2 to 6 to get better results or if code crashes. 
% implement_type -  To detect the gait event, two implementations are provided. Default value: 'fit'
%                   'fit' - 2D Gaussian distribution fitting is done to estimate gait event (runs slow) 
%                   'fast'- faster implemention to estimate the event.
%                   Warning: 'fast' implementation might not give the best
%                   estimate of the gait event compared to gaussian
%                   fitting. It is not presented in the paper but is solely
%                   provided to estimate events faster without gaussian
%                   fitting.

Fs = 128; 
winSizeFactor = 3;
implement_type = 'fit';

% Compute the gait events from the given 3-axis acceleration signal

[HS_LF,TO_LF] = SK_gedAlgo(accX,accY,accZ,Fs,winSizeFactor,implement_type);
                





































    
















