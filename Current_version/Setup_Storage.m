function [spatial] = Setup_Storage(step_cap, C1, C2, C3, myelin)

% SETUP_STORAGE
%
% Preallocates storage for spatial simulation outputs across all timesteps.
% Stores positions of cell populations (C1, C2, C3) and myelin state data.
%
% INPUTS:
%   step_cap - total number of timesteps
%   C1, C2   - initial cell populations
%   C3       - initial (empty or populated) C3 population
%   myelin   - myelin structure containing state information
%
% OUTPUT:
%   spatial  - structure containing preallocated storage arrays

%% Cell Type 1 (C1)
spatial.C1x = cell([1, step_cap + 1]);   % x-positions
spatial.C1y = cell([1, step_cap + 1]);   % y-positions

% Initialise timestep 1
spatial.C1x{1} = C1.x;
spatial.C1y{1} = C1.y;

%% Cell Type 2 (C2)
spatial.C2x = cell([1, step_cap + 1]);   % x-positions
spatial.C2y = cell([1, step_cap + 1]);   % y-positions

% Initialise timestep 1
spatial.C2x{1} = C2.x;
spatial.C2y{1} = C2.y;

%% Cell Type 3 (C3)
spatial.C3x = cell([1, step_cap + 1]);   % x-positions
spatial.C3y = cell([1, step_cap + 1]);   % y-positions

% Initialise timestep 1
spatial.C3x{1} = C3.x;
spatial.C3y{1} = C3.y;

%% Myelin State Storage
spatial.myelinstates = cell([1, step_cap + 1]);  % full state history

% Initialise timestep 1
spatial.myelinstates{1} = myelin.state;

%% Myelin Spatial Categories
spatial.healthy_x  = cell([1, step_cap]);
spatial.healthy_y  = cell([1, step_cap]);

spatial.intermed_x = cell([1, step_cap]);
spatial.intermed_y = cell([1, step_cap]);

spatial.destroy_x  = cell([1, step_cap]);
spatial.destroy_y  = cell([1, step_cap]);

%% Additional State Tracking
spatial.oligo_states  = cell([1, step_cap]);
spatial.myelin_states = cell([1, step_cap]);

end