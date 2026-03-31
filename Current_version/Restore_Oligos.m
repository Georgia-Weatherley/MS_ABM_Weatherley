function myelin = Restore_Oligos( ...
    myelin, replacement_percentage)

% RESTORE_OLIGOS
%
% Reintroduces a proportion of previously non-active (dead) 
% oligodendrocytes and restores their associated myelin to a fully
% healthy state.
%
% INPUTS:
%   myelin                 - myelin structure containing state information
%   replacement_percentage - proportion of dead oligodendrocytes to restore
%
% OUTPUT:
%   myelin - updated myelin structure with restored oligodendrocytes

%% Identify Dead Oligodendrocytes
locate_dead_oligos = find(myelin.oligo_state ~= 1);

%% Determine Number to Restore
num_to_introduce = round( ...
    replacement_percentage * length(locate_dead_oligos));

%% Select Oligodendrocytes to Restore
bring_back = randsample(locate_dead_oligos, num_to_introduce);

% Restore oligodendrocyte functionality
myelin.oligo_state(bring_back) = 1;

%% Restore Associated Myelin
reshaped_myelin = reshape( ...
    myelin.state, ...
    myelin.oligo_dim^2, ...
    myelin.oligo_counter);

% Set restored regions to fully healthy state
reshaped_myelin(:, bring_back) = myelin.grades;

% Reshape back to original structure
myelin.state = reshape(reshaped_myelin, 1, []);

end