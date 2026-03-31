function [xpos, ypos] = Unbiased_Movement(xpos, ypos, PB)

% UNBIASED_MOVEMENT
%
% Performs unbiased random movement of agents on a 2D grid with optional
% probabilistic boundary crossing (e.g. blood-brain barrier).
%
% INPUTS:
%   xpos - x-coordinates of agents
%   ypos - y-coordinates of agents
%   PB   - structure containing boundary crossing probabilities
%
% OUTPUTS:
%   xpos - updated x-coordinates after movement
%   ypos - updated y-coordinates after movement

%% Movement Parameters
delta = 1;                                  % Step size
move_probs = [0.25 0.25 0.25 0.25];         % Equal direction probabilities
cumulative_probs = cumsum(move_probs);

num_cells = length(xpos);

%% Random Sampling
r_C = rand(1, num_cells);   % Random direction selection per agent

%% Westward Movement
west_mask = r_C < cumulative_probs(1);

% Prevent movement beyond left boundary
west_mask(xpos == 1) = 0;

if sum(west_mask) > 0
    r_BBB_west = rand(1, num_cells);
    
    % Acceptance condition for crossing boundary
    accept_mask = r_BBB_west <= PB.enteringprob;
    
    % Restrict crossing if not accepted
    west_mask(west_mask == 1 & xpos == PB.rightside & ...
              accept_mask ~= 1) = 0;
    
    % Apply movement
    xpos(west_mask) = xpos(west_mask) - delta;
end

%% Eastward Movement
east_mask = (r_C < cumulative_probs(2)) & ...
            (r_C >= cumulative_probs(1));

if sum(east_mask) > 0
    r_BBB_east = rand(1, num_cells);
    
    % Acceptance condition for crossing boundary
    accept_mask = r_BBB_east <= PB.leavingprob;
    
    % Restrict crossing if not accepted
    east_mask(east_mask == 1 & xpos == PB.leftside & ...
              accept_mask ~= 1) = 0;
    
    % Apply movement
    xpos(east_mask) = xpos(east_mask) + delta;
end

%% Southward Movement
south_mask = (r_C < cumulative_probs(3)) & ...
             (r_C >= cumulative_probs(2));

ypos(south_mask) = ypos(south_mask) - delta;

%% Northward Movement
north_mask = r_C >= cumulative_probs(3);

ypos(north_mask) = ypos(north_mask) + delta;

end