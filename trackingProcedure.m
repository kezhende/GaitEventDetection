function [sLamda,mu_event,mu_cycle] = trackingProcedure(Fs,rAcc,Wn,scales,nSamp_win,overlap)
% Copyright (c) <2016>, <Siddhartha Khandelwal>

% Tracking procedure to tackle changes in gait speeds
%-------- Input Arguments --------------
%
% Fs    - sampling frequency
% rAcc  - Composite acceleration signal
% Wn    - CWT coefficients
% scales - [smin:smax] the scales chosen for analysis
% nSamp_win - no. of samples in a window
% overlap - overlapping percent
%
%-------- Output Arguments --------------
%
% sLamda - scale that distinguishes event and cycle energies
% mu_event - gaussian mean of the event energies
% mu_cycle - gaussian mean of the cycle energies

plot_tracking = 1; 
if plot_tracking
    h1 = figure; 
end
%% A priori energy density spectrum estimate
x = scales;
mu_er = floor((centfrq('morl')*Fs)./2);
mu_cr = floor((centfrq('morl')*Fs)./1);
a_er = 1; 
a_cr = 1;
sigma_er = round((mu_cr-mu_er)/6); 
sigma_cr = sigma_er;
Es = a_er.*exp(-((x-mu_er)/sigma_er).^2) + a_cr.*exp(-((x-mu_cr)/sigma_cr).^2);

%% Window parameters

if numel(rAcc) <= nSamp_win
    disp('The signal is too short. Enter a longer signal with more samples.')
    return;
end

% Sample nos. for overlapping windows
start_offset = 1:round(nSamp_win*overlap):numel(rAcc)-nSamp_win;

% Preallocate for storing updates
mu_event = zeros(numel(start_offset),1);
mu_cycle = zeros(numel(start_offset),1);
sLamda   = zeros(numel(start_offset),1);
sigma_event = zeros(numel(start_offset),1);
sigma_cycle = zeros(numel(start_offset),1);
a_event = zeros(numel(start_offset),1);
a_cycle = zeros(numel(start_offset),1);
% store_shift = zeros(numel(start_offset),1);

% Init Flag to check if constraint is Not Satisfied
flag_NS = 0;

%% Tracking the Gaussian means/valley in every window

% *************** BEGIN LOOP ********************************************
for i = 1:length(start_offset)
    
    %% Energy density spectrum
    % Select the coefficients
    Wn_r = Wn(:,start_offset(i):start_offset(i)+nSamp_win-1);
    % Compute the energy of the coefficients in the current window
    localEs = abs(Wn_r.*Wn_r); 
    % Normalize
    norm_localEs = localEs./sum(localEs(:));
    % figure; contour(norm_localEs);
    Es_r = sum(norm_localEs,2);
    
    if plot_tracking 
        figure(h1);subplot 511; plot(scales,Es_r,'k'); 
        title('Energy density spectrum'); axis tight;
        subplot 512; plot(scales,Es,'k'); axis tight; 
        title('a priori estimate'); 
    end
        
    %% Cross correlation of a priori estimate with current
    [cc,lag] = xcorr(Es_r,Es);
    [val_cc,index_maxval] = max(cc);
    
    % Measure the scale delay
    tau_r = lag(index_maxval);
   
    if plot_tracking
        % Plot the cross correlation signal
        figure(h1);subplot 513; plot(-(numel(scales)-1):numel(scales)-1,cc,'k'); hold on; 
        plot(tau_r,val_cc,'bo');  xlim([-150 150]);
        title('$E_{s,r}^- \star E_{s,r}$','Interpreter','latex'); 
        set(gca,'xgrid','on'); hold off;
    end
    
   %% CONSTRAINT CHECK
    
   % CONSTRAINT 1: 
   if mu_er+tau_r > 1 && mu_cr+tau_r < scales(end)
        % compute scale that separates event and cycle energies
        [~,idx_sLamda_r] = min(Es_r(mu_er+tau_r:mu_cr+tau_r));
        sLamda_r = (mu_er+tau_r) + idx_sLamda_r - 1;
        
        % Compute the aposteriori estimate parameters
        [~,mu_er] = max(Es_r(1:sLamda_r)); 
        [~,idx_mu_cr] = max(Es_r(sLamda_r:end)); 
        mu_cr = idx_mu_cr + sLamda_r - 1;        
        a_er = Es_r(mu_er); 
        a_cr = Es_r(mu_cr);
        
        % A posteriori estimate
        Es =  a_er.*exp(-((x-mu_er)/sigma_er).^2) + a_cr.*exp(-((x-mu_cr)/sigma_cr).^2);
   else
        % Set Flag if the constraint is violated
        flag_NS = 1;
   end
     
   % CONSTRAINT 2:
   if (mu_cr/mu_er) < 1.9 || (mu_cr/mu_er) > 2.1       
       % Set Flag if the constraint is violated
       flag_NS = 1;
   end
   
    %% Constraints not satisfied - fit a two term 1D gauss mix
    
    if flag_NS == 1
        [a_er,mu_er,sigma_er,a_cr,mu_cr,sigma_cr] = bestFit(x,Es_r);
       
        % Relax Constraint II more to accommodate fit. If violated, copy parameters
         if ((mu_cr/mu_er) < 1.8 || (mu_cr/mu_er) > 2.2) && i > 1
             a_er = a_event(i-1);
             a_cr = a_cycle(i-1);
             mu_er = mu_event(i-1);
             mu_cr = mu_cycle(i-1);
             sigma_er = sigma_event(i-1);
             sigma_cr = sigma_cycle(i-1);
             sLamda_r = sLamda(i-1);
             Es = a_er.*exp(-((x-mu_er)/sigma_er).^2) + a_cr.*exp(-((x-mu_cr)/sigma_cr).^2);
         else
             %  Compute sLamda and a posteriori estimate
             [~,sLamda_r] = min(Es_r(mu_er:mu_cr));
             sLamda_r = mu_er + sLamda_r - 1;
             Es = a_er.*exp(-((x-mu_er)/sigma_er).^2) + a_cr.*exp(-((x-mu_cr)/sigma_cr).^2);
         end
        
        % Reset the Flag to 0 
        flag_NS = 0;
    end
    
    %% Store the updated gaussian means 
    mu_event(i) = mu_er;
    mu_cycle(i) = mu_cr;
    sLamda(i) = sLamda_r;
    sigma_event(i) = sigma_er;
    sigma_cycle(i) = sigma_cr;
    a_event(i) = a_er;
    a_cycle(i) = a_cr;    
    
    if plot_tracking
        % Plot the updated gaussian function
        figure(h1);
        subplot 514; plot(scales,Es,'k'); hold on; 
        title('\emph{a posteriori} estimate $\hat{E}_{s,r}$','Interpreter','latex'); 
        axis tight; hold off;
        % Plot the identified valley and new gaussain means
        subplot 515; plot(Es_r); hold on;
        title('Tracking the Velocity'); axis tight;  xlabel('Scales');
        plot(sLamda(i),Es_r(sLamda(i)),'bo');
        plot([mu_event(i),mu_cycle(i)],[Es_r(mu_event(i)),Es_r(mu_cycle(i))],'ro'); hold off;
        pause(0.001)
    end    
end
% *************** END LOOP ********************************************

end



