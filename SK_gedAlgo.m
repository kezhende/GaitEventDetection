function [vec_HS,vec_TO] = SK_gedAlgo(accX,accY,accZ,Fs,winSizeFactor,implement_type)
% Copyright (c) <2016>, <Siddhartha Khandelwal>

% Publication to cite: 
% Siddhartha Khandelwal, Nicholas Wickström (2016). 
% Gait Event Detection in Real-World Environment for Long-Term Applications : 
% Incorporating Domain Knowledge into Time-Frequency Analysis. 
% IEEE transactions on neural systems and rehabilitation engineering.
%-------- Input Arguments --------------
%
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
% 
%-------- Output Arguments --------------
%
% HS - vector containing sample no. of Heel Strike events
% TO - vector containing sample no. of Toe Off events
%% Input checks
if nargin < 5
    error('Not enough input arguments');
else if nargin == 5
        implement_type = 'fit';
    end
end
%% Composite acceleration signal
rAcc = sqrt(accX.^2 + accY.^2 + accZ.^2);
% figure; plot(rAcc);

%% CWT 

scales = (1:floor((centfrq('morl')*Fs)./0.5))';
Wn = cwt(rAcc,scales,'morl');

%% Tracking Procedure
% No. of samples in a window
nSamp_win = floor(winSizeFactor*Fs); 
overlap = (1 - 50/100);
[sLamda,mu_e,mu_c] = trackingProcedure(Fs,rAcc,Wn,scales,nSamp_win,overlap);

%% Obtaining temporal signals 

[x_en,x_cn,vec_mu_e,~,sLamda_n] = temporalSignal(mu_e,mu_c,sLamda,rAcc,Wn,nSamp_win,overlap,Fs);
% figure; plot(zscore(rAcc)); hold on; plot(zscore(x_en),'r'); plot(zscore(x_cn),'k');

%% Locating 2-D spectral temporal event regions 

% Local minima points in x_en
[~,m_local] = findpeaks(-x_en);
% figure; plot(zscore(x_en)); hold on; plot(m_local,-ones(numel(m_local)),'go');

% Add the first and last sample of the signal if necessary
if m_local(1) > median(diff(m_local)) - 10
    m_local = [1; m_local];
else if numel(rAcc) - m_local(end) > median(diff(m_local)) - 10
        m_local = [m_local; numel(rAcc)];
    end
end

%% Detect the gait event

% Event Detection: 2D Gaussian Fitting
gaitEvent = detectEvent(Wn,sLamda_n,m_local,vec_mu_e,implement_type);

% figure; plot(zscore(sig)); hold on; plot(gaitEvent,1.5.*ones(size(gaitEvent)),'ko');

%% Label the event as HS or TO

% Event Detection: 2D Gaussian Fitting
[vec_HS,vec_TO] = labelEvent(x_cn,gaitEvent,rAcc);


%% Final Plots
% % Plot the stride time
% figure;subplot(2,1,1);plot(diff(vec_HS)./Fs,'k'); title('HS');
% subplot(2,1,2);plot(diff(vec_TO)./Fs,'k');title('TO');


% % % Plot the signal with the detected HS and TO
figure;
plot(rAcc); hold on;
plot(vec_HS,40*ones(length(vec_HS),1),'go');
plot(vec_TO,40*ones(length(vec_TO),1),'k^');
legend('rAcc Signal','HS','TO');
xlabel('Time [samples]');
ylabel('Amplitude [m/s^2]');
end

