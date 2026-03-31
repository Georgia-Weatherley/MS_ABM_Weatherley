function [spatial] = Tidy_To_Plot(spatial)

% TIDY_TO_PLOT
%
% Fills empty cells in spatial data structures with Inf for plotting purposes.
%
% INPUT:
%   spatial - structure containing cell arrays for C1, C2, C3, and myelin states
%
% OUTPUT:
%   spatial - updated structure with no empty cells

%% Replace empty fields with inf 

spatial.C1x(cellfun('isempty',spatial.C1x)) = {Inf}; 
spatial.C1y(cellfun('isempty',spatial.C1y)) = {Inf};
spatial.C2x(cellfun('isempty',spatial.C2x)) = {Inf}; 
spatial.C2y(cellfun('isempty',spatial.C2y)) = {Inf};
spatial.C3x(cellfun('isempty',spatial.C3x)) = {Inf}; 
spatial.C3y(cellfun('isempty',spatial.C3y)) = {Inf};
spatial.healthy_x(cellfun('isempty',spatial.healthy_x)) = {Inf}; 
spatial.healthy_y(cellfun('isempty',spatial.healthy_y)) = {Inf};
spatial.intermed_x(cellfun('isempty',spatial.intermed_x)) = {Inf}; 
spatial.intermed_y(cellfun('isempty',spatial.intermed_y)) = {Inf};
spatial.destroy_x(cellfun('isempty',spatial.destroy_x)) = {Inf}; 
spatial.destroy_y(cellfun('isempty', spatial.destroy_y)) = {Inf};

end