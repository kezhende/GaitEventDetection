function [vec_HS,vec_TO] = labelEvent(x_cn,gaitEvent,rAcc)
% Label the detected gait events as HS or TO
%  
%-------- Input Arguments --------------
%
% x_cn - cycle temporal signal
% gaitEvent - sample nos. of the estimated gait events 
% rAcc - Composite acceleration signal
%
%-------- Output Arguments --------------
%
% vec_HS - vector containing sample no. of Heel Strike events
% vec_TO - vector containing sample no. of Toe Off events
%% 

% Make a copy of the cycle cross section curve
copy_x_cn = x_cn;
% Set all positive values to one and negative to zero
copy_x_cn(copy_x_cn > 0) = 1;
copy_x_cn(copy_x_cn < 0) = 0;
% The cycle starts when the curve goes from negative to positive 
start_cycle = find(diff(copy_x_cn) == 1);
% Define the cycle boundaries
cycle_bound(:,1) = [1;start_cycle];
cycle_bound(:,2) = [start_cycle;numel(rAcc)];

% Initialize the HS and TO vectors
vec_HS = zeros(size(cycle_bound,1),1);
vec_TO = zeros(size(cycle_bound,1),1);

for i = 1:size(cycle_bound,1)
    
    % Events that lie in the current cycle
    x1 = find(gaitEvent >= cycle_bound(i,1));
    x2 = find(gaitEvent <= cycle_bound(i,2));
    % Index of the events that lie in the current cycle
    index_currentevent = intersect(x1,x2);
    
    %  Generally two events should lie within a cycle
    if numel(index_currentevent) == 2
        % First event is labelled as HS
        vec_HS(i) = gaitEvent(index_currentevent(1));
        % Second event is labelled as TO
        vec_TO(i) = gaitEvent(index_currentevent(2));
    end
    
    % If no. of event in a cycle ~= 2
    if numel(index_currentevent) == 1
        % Assume event as a HS
        vec_HS(i) = gaitEvent(index_currentevent);
    else if numel(index_currentevent) > 2
            % Consider first and second events in the cycle
            vec_HS(i) = gaitEvent(index_currentevent(1));
            vec_TO(i) = gaitEvent(index_currentevent(2));
        end
    end    
end

% Remove any zeros
vec_HS (vec_HS == 0) = [];
vec_TO (vec_TO == 0) = [];

%% Plots
% figure; plot(zscore(rAcc)); hold on; plot(zscore(x_cn));
% plot(start_cycle,zeros(length(start_cycle)),'ko');
% plot(vec_HS,4*ones(length(vec_HS),1),'go');
% plot(vec_TO,4*ones(length(vec_TO),1),'m^');

end

