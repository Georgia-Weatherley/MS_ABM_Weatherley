function [result, tracked, myelin, spatial] = main_function( ...
    save_spatial, intervention, relapse_schedule, step_cap)

% MAIN_FUNCTION
%
% Runs a single stochastic simulation of the agent-based model.
%
% INPUTS:
%   save_spatial     - (logical) check if saving spatial data
%   intervention     - intervention parameters and switches
%   relapse_schedule - relapse schedule identifier
%   step_cap - number of time steps for a single realisation 
%
% OUTPUTS:
%   result  - population-level time series data
%   tracked - tracked myelin and oligodendrocyte states
%   myelin  - final myelin state
%   spatial - spatial outputs (if enabled)

%% Initial Parameters
threshold_oligo_apop      = intervention.oligo_apop;
threshold_stop_myelination = intervention.oligo_stop_my;

%% Setup Domain and Initial Conditions 
[funcs, PB, sizing, borders, myelin] = Setup_Simulation;
[myelin, C1, C2, C3, tau] = Setup_Agents(myelin, borders, sizing);

%% Tracking Setup
% Predefined indices for tracking specific myelin/oligos

indx1 = [977 935 1036]; % oligo indices
indx2 = [133 168 232]; % oligo indices
indx3 = [429 486 529]; % oligo indices

m_indx1 = [24413 23363 25888]; % myelin indices
m_indx2 = [3313 4188 5788];  % myelin indices
m_indx3 = [10713 12138 13213];  % myelin indices

tracked.indices        = [m_indx1 m_indx2 m_indx3]; 
tracked.oligo_indices  = [indx1 indx2 indx3];

%% Storage Allocation
tracked.states        = zeros(length(tracked.indices), step_cap);
tracked.oligo_states  = zeros(length(tracked.oligo_indices), step_cap);

result.myelinating_oligo     = zeros(1, step_cap);
result.dead_oligo            = zeros(1, step_cap);
result.nonmyelinating_oligo  = zeros(1, step_cap);

result.healthy_myelin        = zeros(1, step_cap);
result.intermediate_myelin   = zeros(1, step_cap);
result.destroyed_myelin      = zeros(1, step_cap);

result.C1_pop = zeros(1, step_cap);
result.C2_pop = zeros(1, step_cap);
result.C3_pop = zeros(1, step_cap);

if save_spatial == true
    spatial = Setup_Storage(step_cap, C1, C2, C3, myelin);
else
    spatial = 0;
end

%% Main Time Loop

for nn = 1:step_cap
    
    %% Apply Interventions
    if intervention.BBB == "on" && (sum(nn == intervention.BBB_times) == 1) % time to change BBB permeability 
        PB.leavingprob = intervention.new_BBB_prob; % update permeability 
    end

    if intervention.restore_oligos == "on" && ...
            (sum(nn == intervention.new_oligo_times) == 1) % time to restore oligos 
        myelin = Restore_Oligos( ...
            myelin, intervention.restore_proportion); % restore oligos 
    end

    if intervention.oligo_properties == "on" && ...
            (sum(nn == intervention.oligo_prop_times) == 1) % time to change oligo properties 
        threshold_oligo_apop       = intervention.new_apop; % change apoptosis resilience
        threshold_stop_myelination = intervention.new_stopmy; % change stop remyelination resilience
    end
    
    %% Cell Movement
    % C3 movement (biased/unbiased depending on myelin state)
    if length(C3.x) ~= 0
        myelin.sum = sum(myelin.state);
        
        if myelin.sum == myelin.fullyhealthyquant || myelin.sum == 0 % no myelin left or all myelin intact 
            [C3.x, C3.y] = Unbiased_Movement(C3.x, C3.y, PB); % unbiased movement of C3 
        else
            C3 = Biased_Movement(C3, myelin, funcs, PB); % biased movement of C3
        end
    end

    % C1 and C2 movement (unbiased)
    [C1.x, C1.y] = Unbiased_Movement(C1.x, C1.y, PB);
    [C2.x, C2.y] = Unbiased_Movement(C2.x, C2.y, PB);

    %% Boundary Conditions
    % Vertical wrapping
    wrap_down.C1 = C1.y == sizing.domainheight; C1.y(wrap_down.C1) = 1;
    wrap_down.C2 = C2.y == sizing.domainheight; C2.y(wrap_down.C2) = 1;
    wrap_down.C3 = C3.y == sizing.domainheight; C3.y(wrap_down.C3) = 1;

    wrap_up.C1 = C1.y == 0; C1.y(wrap_up.C1) = sizing.domainheight - 1;
    wrap_up.C2 = C2.y == 0; C2.y(wrap_up.C2) = sizing.domainheight - 1;
    wrap_up.C3 = C3.y == 0; C3.y(wrap_up.C3) = sizing.domainheight - 1;

    % Horizontal boundaries
    C1.x(C1.x == 0) = 1;  % Reflecting boundary (LHS)
    C2.x(C2.x == borders.rightcns + 1) = borders.rightcns;  % RHS

    % Eviction masks
    evict.C1 = C1.x == sizing.domainwidth | C1.x == 0;
    evict.C3 = C3.x == sizing.domainwidth;

    % Remove evicted cells
    C1.x = C1.x(~evict.C1); C1.y = C1.y(~evict.C1);
    C3.x = C3.x(~evict.C3); C3.y = C3.y(~evict.C3);

    %% Myelin Dynamics
    myelin = Oligo_Remyelination( ...
        threshold_stop_myelination, myelin); % remyelinate 

    myelin = Oligo_Apoptosis( ...
        myelin, threshold_oligo_apop, threshold_stop_myelination); % check which oligodendrocytes exceed thresholds 

    if length(C3.x) ~= 0
        myelin = Demyelination(myelin, C3); % check for demyelination by C3 cells 
    end

    %% Cell Interactions (C1 + C2 create C3)
    if length(C1.x) ~= 0 & length(C2.x) ~= 0
        C1.pos = [C1.x; C1.y]'; % get positions of C1 cells 
        C2.pos = [C2.x; C2.y]'; % get positions of C2 cells 

        [~, index_C1, ~] = intersect(C1.pos, C2.pos, 'rows'); % find where they occupy the same lattice site 

        new_C3x_pos = C1.x(index_C1); % add C3 cell to lattice site 

        C3.x = [C3.x new_C3x_pos]; % update array tracking C3 positions 
        C3.y = [C3.y C1.y(index_C1)];

        C1.x = C1.x(setdiff(1:end, index_C1)); % remove C1 from tracker to acknowledge conversion to C3
        C1.y = C1.y(setdiff(1:end, index_C1));
    end

    %% Add C1, Cell Death
    C1 = Renew_Cell_1(C1, nn, tau, relapse_schedule); % add new C1 cells according to relapse schedule 

    mort.C1mask = rand(1, length(C1.x)) < C1.deathrate; 
    mort.C3mask = rand(1, length(C3.x)) < C3.deathrate;

    C1.x = C1.x(~mort.C1mask); C1.y = C1.y(~mort.C1mask); % remove C1 cells due to apoptosis 
    C3.x = C3.x(~mort.C3mask); C3.y = C3.y(~mort.C3mask); % remove C3 cells due to apoptosis 

    %% Data Collection
    result.myelinating_oligo(nn)    = sum(myelin.oligo_state == 1);
    result.dead_oligo(nn)           = sum(myelin.oligo_state == 0);
    result.nonmyelinating_oligo(nn) = sum(myelin.oligo_state == 2);

    result.healthy_myelin(nn)      = sum(myelin.state == myelin.grades);
    result.intermediate_myelin(nn) = sum(myelin.state > 0 & ...
                                         myelin.state < myelin.grades);
    result.destroyed_myelin(nn)    = sum(myelin.state == 0);

    result.C1_pop(nn) = length(C1.x);
    result.C2_pop(nn) = length(C2.x);
    result.C3_pop(nn) = length(C3.x);

    tracked.states(:, nn)        = myelin.state(tracked.indices);
    tracked.oligo_states(:, nn)  = ...
        myelin.oligo_state(tracked.oligo_indices);

    %% Spatial Data Storage
    if save_spatial == true
        
        spatial.C1x{nn+1} = C1.x; spatial.C1y{nn+1} = C1.y;
        spatial.C2x{nn+1} = C2.x; spatial.C2y{nn+1} = C2.y;
        spatial.C3x{nn+1} = C3.x; spatial.C3y{nn+1} = C3.y;

        spatial.healthy_x{nn}  = ...
            myelin.xarray(myelin.state == myelin.grades);
        spatial.healthy_y{nn}  = ...
            myelin.yarray(myelin.state == myelin.grades);

        spatial.intermed_x{nn} = ...
            myelin.xarray(myelin.state < myelin.grades & myelin.state > 0);
        spatial.intermed_y{nn} = ...
            myelin.yarray(myelin.state < myelin.grades & myelin.state > 0);

        spatial.destroy_x{nn}  = ...
            myelin.xarray(myelin.state == 0);
        spatial.destroy_y{nn}  = ...
            myelin.yarray(myelin.state == 0);

        spatial.oligo_states{nn}  = myelin.oligo_state;
        spatial.myelin_states{nn} = myelin.state;
    end
end

%% Post-processing
if save_spatial == true
    spatial = Tidy_To_Plot(spatial);
end
