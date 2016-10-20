function [filt_x_en,filt_x_cn,vec_mu_e,vec_mu_c,sLamda_n] = temporalSignal(mu_e,mu_c,sLamda,rAcc,Wn,nSamp_win,overlap,Fs)
% Copyright (c) <2016>, <Siddhartha Khandelwal>

% Obtain temporal signals that match the frequency of event and cycle 
%-------- Input Arguments --------------
%
% mu_e,mu_c,sLamda - parameters from tracking procedure
% rAcc - Composite acceleration signal
% Wn - CWT coefficients
% nSamp_win - no. of samples in a window
% overlap - overlapping percent
%
%-------- Output Arguments --------------
% filt_x_en,filt_x_cn - low-pass filtered event and cycle temporal signals 
% vec_mu_e,vec_mu_c - gaussian means of event and cycle energies, compiled
%                       from all windows
% sLamda_n - sLamda complied from all windows
%%

% Resize the vectors from tracking procedure
vec_mu_e = repmat(mu_e',round(nSamp_win*overlap),1);
vec_mu_c = repmat(mu_c',round(nSamp_win*overlap),1);
sLamda_n = repmat(sLamda',round(nSamp_win*overlap),1);

vec_mu_e = reshape(vec_mu_e,size(vec_mu_e,1)*size(vec_mu_e,2),1);
vec_mu_c = reshape(vec_mu_c,size(vec_mu_c,1)*size(vec_mu_c,2),1);
sLamda_n = reshape(sLamda_n,size(sLamda_n,1)*size(sLamda_n,2),1);

if numel(vec_mu_e) < numel(rAcc)
    vec_mu_e = [vec_mu_e; vec_mu_e(end).*ones(numel(rAcc)-numel(vec_mu_e),1)];
    vec_mu_c = [vec_mu_c; vec_mu_c(end).*ones(numel(rAcc)-numel(vec_mu_c),1)];
    sLamda_n = [sLamda_n; sLamda_n(end).*ones(numel(rAcc)-numel(sLamda_n),1)];
end

% Obtain the temporal signals
idx_event = sub2ind(size(Wn), vec_mu_e', 1:numel(vec_mu_e));
idx_cycle = sub2ind(size(Wn), vec_mu_c', 1:numel(vec_mu_c));
x_en = Wn(idx_event)';
x_cn = Wn(idx_cycle)';

% Remove high freq. noise from the temporal signals
% Low pass FIR filter with Cut-off freq = 8Hz
cutoffFreq = 8;
lpFilt = designfilt('lowpassfir','PassbandFrequency',cutoffFreq/(Fs/2), ...
         'StopbandFrequency',(cutoffFreq+1)/(Fs/2),'PassbandRipple',0.5, ...
         'StopbandAttenuation',60,'DesignMethod','equiripple');
%fvtool(lpFilt)
filt_x_en = filtfilt(lpFilt,x_en);
filt_x_cn = filtfilt(lpFilt,x_cn);


%% Plots

% figure; plot(zscore(rAcc)); hold on; plot(zscore(filt_x_en),'r'); plot(zscore(filt_x_cn),'k');

end

