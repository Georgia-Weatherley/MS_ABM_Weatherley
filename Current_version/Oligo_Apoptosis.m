function myelin = Oligo_Apoptosis(myelin, threshold_apop, threshold_stopmy)

% OLIGO_APOPTOSIS
%
% Implements oligodendrocyte apoptosis and function loss based on
% accumulated damage in myelin segments.
%
% INPUTS:
%   myelin            - structure containing myelin state and oligodendrocyte info
%   threshold_apop    - damage threshold for apoptosis (cell death)
%   threshold_stopmy  - damage threshold for halting remyelination
%
% OUTPUT:
%   myelin - updated structure with modified oligodendrocyte states
%
% NOTE:
%   - Oligos with damage >= threshold_stopmy stop myelinating (state = 2)
%   - Oligos with damage >= threshold_apop undergo apoptosis (state = 0)

%% Reshape Myelin by Oligodendrocyte
reshaped_myelin = reshape(myelin.state, ...
                          myelin.oligo_dim^2, ...
                          myelin.oligo_counter);

%% Compute Damage per Oligo
indicate_destroyed_myelin = reshaped_myelin == 0;
total_damage_to_oligo = sum(indicate_destroyed_myelin);

%% Identify Oligos That Stop Myelinating
index_stopmy = find(total_damage_to_oligo >= threshold_stopmy);
myelin.oligo_state(index_stopmy) = 2; % no longer able to remyelinate

%% Identify Oligos That Undergo Apoptosis
index_apop = find(total_damage_to_oligo >= threshold_apop);
myelin.oligo_state(index_apop) = 0;       % oligodendrocyte is dead
reshaped_myelin(:, index_apop) = 0;       % remove all associated myelin

%% Reshape Myelin State
myelin.state = reshape(reshaped_myelin, 1, []);

end