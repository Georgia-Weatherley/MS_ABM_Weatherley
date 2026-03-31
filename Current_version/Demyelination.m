function [myelin] = Demyelination(myelin, C3)

% DEMYELINATION
%
% Degrades myelin segments that are in contact with C3 cells.
%
% INPUTS:
%   myelin - structure containing myelin state, positions, and timers
%   C3     - structure containing C3 cell positions (x and y arrays)
%
% OUTPUT:
%   myelin - updated myelin structure with degraded states
%
% NOTES:
%   - Myelin state is reduced by 1 grade where occupied by C3
%   - Timer for degraded myelin is reset to -1

%% Prepare C3 Positions
C3_pos = [C3.x; C3.y]';

%% Find Myelin Occupied by C3
[~, index_myelin, ~] = intersect(myelin.colarray, C3_pos, 'rows');

%% Degrade Myelin State
myelin.state(index_myelin) = myelin.state(index_myelin) - 1;

%% Reset Timer for Degraded Myelin
myelin.timer(index_myelin) = -1;

%% Ensure No Negative Myelin States (in case of simulataneous degradation)
myelin.state(myelin.state < 0) = 0;