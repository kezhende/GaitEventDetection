function [a_e,mu_e,sigma_e,a_c,mu_c,sigma_c] = bestFit(x,Es)
% Copyright (c) <2016>, <Siddhartha Khandelwal>

% Two-term 1D gaussian fitting over Es
%-------- Input Arguments --------------
% Input
% x - Scales
% Es - Energy density spectrum estimate

%-------- Output Arguments --------------
% a_e, mu_e, sigma_e - gauss parameters: event 
% a_c, mu_c, sigma_c - gauss parameters: cycle 


%%
% Most energetic scale in Es
[a_x,mu_x] = max(Es);

% Fit Options
Set1 = [a_x, mu_x/2, x(end)/16, a_x, mu_x,   x(end)/8];
Set2 = [a_x, mu_x,   x(end)/16, a_x, 2*mu_x, x(end)/8];       
% First Set
fo1 = fitoptions('gauss2');
fo1.Robust = 'Bisquare';
fo1.StartPoint = Set1;
% Second Set
fo2 = fitoptions('gauss2');
fo2.Robust = 'Bisquare';
fo2.StartPoint = Set2;

% Fit atwo term 1D gaussian
gaussfit1 = fit(x,Es,'gauss2',fo1);
gaussfit2 = fit(x,Es,'gauss2',fo2);

%% Compute RMSE for both fits

y1 = feval(gaussfit1,x);
y2 = feval(gaussfit2,x);
idx_y1 = find(y1 > max(Es)/10);
idx_y2 = find(y2 > max(Es)/10);

rmse1 = sqrt(mean((y1(idx_y1) - Es(idx_y1)).^2));
rmse2 = sqrt(mean((y2(idx_y2) - Es(idx_y2)).^2));

%% Evaluate the best fit

% Init Flag to check if bounds are Not Satisfied
flag_NS = [0,0];
% Check if the fit parameters lie within the expected bounds
fitParameters = [coeffvalues(gaussfit1) ; coeffvalues(gaussfit2)];
for i = 1:2
    a = [fitParameters(i,1) fitParameters(i,4)];
    mu = [fitParameters(i,2) fitParameters(i,5)];
    sigma = [fitParameters(i,3) fitParameters(i,6)];
    % amplitude bounds
    if ~isempty(a(a < 0))
        flag_NS(i) = 1;
        % mean bounds
    else if ~isempty(mu(mu < 0) | mu(mu > x(end)))
            flag_NS(i) = 1;
            % sigma bounds
        else if ~isempty(sigma(sigma < 0) | sigma(sigma > x(end)))
                flag_NS(i) = 1;
            end
        end
    end
end
            
if ~isempty(flag_NS(flag_NS == 1))
    % Discard the fit if a flag was set for a fit
    if flag_NS(1) == 1       
        bFitParams = fitParameters(2,:);
    else
        bFitParams = fitParameters(1,:);
    end
else
    % If no flags were set, compare rmse of both fits
    if rmse1 < rmse2
        bFitParams = fitParameters(1,:);
    else
        bFitParams = fitParameters(2,:);
    end
end

%% Assign the respective gaussian parameters based on the mu

if bFitParams(2) < bFitParams(5)    
    a_e = bFitParams(1);
    mu_e = floor(bFitParams(2));
    sigma_e = bFitParams(3);
    a_c = bFitParams(4);
    mu_c = floor(bFitParams(5));
    sigma_c = bFitParams(6);
else
    a_c = bFitParams(1);
    mu_c = floor(bFitParams(2));
    sigma_c = bFitParams(3);
    a_e = bFitParams(4);
    mu_e = floor(bFitParams(5));
    sigma_e = bFitParams(6);
end

%% Plot 

% figure; plot(Es); hold on; plot(gaussfit1,'r'); plot(gaussfit2,'k'); hold on;
% plot(idx_y1,y1(idx_y1),'or'); plot(idx_y2,y2(idx_y2),'ok');
% plot(a_e.*exp(-((x-mu_e)/sigma_e).^2) + a_c.*exp(-((x-mu_c)/sigma_c).^2),'g','Linewidth',1.5);
end
    








