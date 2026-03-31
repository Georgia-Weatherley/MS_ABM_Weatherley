function [myelin] = Oligo_Remyelination(local_threshold, myelin)

% OLIGO_REMYELINATION
%
% Simulates remyelination driven by oligodendrocytes.
% Myelin regeneration occurs if:
%   - Local damage is below a threshold
%   - The associated oligodendrocyte is healthy
%   - A regeneration timer condition is met
%
% INPUTS:
%   local_threshold - damage threshold for regeneration
%   myelin          - structure containing myelin state and metadata
%   funcs           - structure containing regeneration functions (unused here)
%
% OUTPUT:
%   myelin - updated structure after remyelination step

%% Reshape Myelin State
reshaped_myelin = reshape(myelin.state, ...
                          myelin.oligo_dim^2, ...
                          myelin.oligo_counter);

%% Compute Damage per Oligo
indicate_destroyed_myelin = reshaped_myelin == 0;

% Total damage per oligodendrocyte (column-wise)
total_damage_to_oligo = sum(indicate_destroyed_myelin);

% Expand to match full grid
dim_fix_damage = repmat(total_damage_to_oligo, ...
                        myelin.oligo_dim^2, 1);

%% Identify Healthy Oligodendrocytes
oligo_healthy = myelin.oligo_state == 1;

dim_fix_oligo = repmat(oligo_healthy, ...
                       myelin.oligo_dim^2, 1);

array_healthy_oligos = reshape(dim_fix_oligo, 1, []);

%% Update Regeneration Timers
% Increase timer for myelin associated with healthy oligodendrocytes
myelin.timer(array_healthy_oligos) = ...
    myelin.timer(array_healthy_oligos) + 1;

%% Regeneration Condition
timer_mask = myelin.timer >= myelin.healtime;

dim_fix_timer = reshape(timer_mask, ...
                        myelin.oligo_dim^2, ...
                        myelin.oligo_counter);

%% Identify Myelin to Regenerate
build_back_mask = ...
    (dim_fix_damage < local_threshold) & ...
    (dim_fix_timer == 1) & ...
    (dim_fix_oligo == 1);

array_build_back_mask = reshape(build_back_mask, 1, []);

%% Apply Remyelination
myelin.state(array_build_back_mask) = ...
    myelin.state(array_build_back_mask) + 1;

% Reset timers for regenerated myelin
myelin.timer(array_build_back_mask) = 0;

%% Enforce Bounds
% Prevent exceeding maximum health
myelin.state(myelin.state > myelin.grades) = myelin.grades;

% Fully healthy myelin no longer regenerates
myelin.timer(myelin.state == myelin.grades) = -1;

end