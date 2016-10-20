function [gaitEvent] = detectEvent(Wn,sLamda_n,m_local,vec_mu_e,implement_type)
% Copyright (c) <2016>, <Siddhartha Khandelwal>

% Estimate the gait event from 2D spectral temporal event region
%
%-------- Input Arguments --------------
%
% Wn - CWT coefficients
% sLamda_n - sLamda complied from all windows
% m_local - Local minima points in x_en
% vec_mu_e - gaussian means of event energies compiled from all windows
% implement_type - 'fit': 2D gaussian surface fitting over event region
%                - 'fast': a faster algo to estimate event without fitting

%-------- Output Arguments --------------
% gaitEvent - sample nos. of the estimated gait events 

%%
% Initialize vector
gaitEvent = zeros(length(m_local)-1,1);

if strcmp(implement_type,'fit')

    % General model for least-squares surface gaussian distribution fit - Streutker 2003    
    % This general model is not a part of the paper but it gives better results
    % ft(a0,a_e,mu_x,mu_y,phi,sigma_x,sigma_y,x,y) 
    ft = fittype( 'a0 + a_e*exp(-(((x-mu_x)*cosd(phi)+(y-mu_y)*sind(phi))^2)/(0.5*(sigma_x^2))-(((y-mu_y)*cosd(phi)+(x-mu_x)*sind(phi))^2)/(0.5*(sigma_y^2)))','independent',{'x', 'y'},'dependent','z'); 
    fop = fitoptions(ft);
    fop.Lower = [0; 0; 0 ; 0; 0; 0; 0];

    for i = 1:numel(m_local)-1   
        % Select the sLamda bound
        sLamda_bound = round(mean(sLamda_n(m_local(i):m_local(i+1))));
        % 2D spectral temporal event region
        R_ns = Wn(1:sLamda_bound,m_local(i):m_local(i+1));    

        % Declare fitoptions : StartPoint and Upper.
        mu_e = round(mean(vec_mu_e(m_local(i):m_local(i+1))));    
        fop.Upper = [50; max(max(R_ns)); size(R_ns,2); size(R_ns,1); 90; size(R_ns,2)/2; size(R_ns,1)/2];
        fop.StartPoint = [0; max(max(R_ns)); size(R_ns,2)/2; mu_e; 0; 5; 5];
        % Prepare coeffs for surface fitting
        [xGrid,yGrid]=meshgrid(1:size(R_ns,2),1:size(R_ns,1));
        [xOut, yOut, zOut] = prepareSurfaceData(xGrid,yGrid,R_ns);
        % Fit   
        fitresult = fit( [xOut, yOut], zOut, ft, fop );

        % Store where the event happened
        gaitEvent(i) = round(fitresult.mu_x) + m_local(i) - 1;          
    end
    
else if strcmp(implement_type,'fast')
                
        for i = 1:numel(m_local)-1   
            % Select the sLamda bound
            sLamda_bound = round(mean(sLamda_n(m_local(i):m_local(i+1))));
            % 2D spectral temporal event region
            R_ns = Wn(1:sLamda_bound,m_local(i):m_local(i+1));  
            % Maximum value along the frequency and time axis      
            [~,cc_freq] = max(R_ns);    
            [~,cc_time] = max(R_ns,[],2); 
            
            %Find the intersection of cc_time and cc_freq curves
            dmatrix = 99.*ones(size(R_ns));
            browse = sort(unique(cc_time))';
            for k = 1:length(browse)
                rt = find(cc_time == browse(k));
                idx = sub2ind(size(dmatrix),rt, browse(k).*ones(size(rt)));
                dmatrix(idx) = abs(rt - cc_freq(browse(k)));
            end
            % Minimum value in dmatrix i.e. 0 & 1 is where the intersection is
            common = ismember(dmatrix,[0 1]);
            [intersectFreq_pts,intersectTime_pts] = find(common);
            if length(intersectFreq_pts) > 1 
                % Select the mu_event for the cuurent window 
                mu_e = round(mean(vec_mu_e(m_local(i):m_local(i+1))));                 
                % Find the freq_pt that is closest to mu_event
                [~,idx_close] = min(abs(intersectFreq_pts - mu_e));
                % Store the intersection point as gait event in time
                gaitEvent(i) = intersectTime_pts(idx_close) + m_local(i) - 1;                    
            else if length(intersectFreq_pts) == 1 
                    % If there is only one intersection
                    gaitEvent(i) = intersectTime_pts + m_local(i) - 1; 
                end
            end            
        end        
    end
end




