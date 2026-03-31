function [C1] = Renew_Cell_1(C1, nn, tau, relapse_schedule)

% RENEW_CELL_1
%
% Introduces new C1 cells (T cells) in the peripheral blood according to
% specified relapse schedules.
%
% INPUTS:
%   C1              - structure containing current C1 cells (x, y, age)
%   nn              - current timestep
%   tau             - timestep duration (minutes)
%   relapse_schedule- integer (1-4) specifying relapse schedule type
%
% OUTPUT:
%   C1 - updated structure with new cells added
%
% NOTES:
%   - Supports four relapse schedules (3 or 6 relapses of varying duration)

%% Setup
mins_per_day = 1440;
steps_per_day = mins_per_day / tau;

%% Define Relapse Windows
switch relapse_schedule
    case 1 % 3 relapses, 28 days each
        days_relapse = 28;
        start_days = [1, 100, 200];
    case 2 % 3 relapses, 60 days each
        days_relapse = 60;
        start_days = [1, 100, 200];
    case 3 % 6 relapses, 14 days each
        days_relapse = 14;
        start_days = [1, 50, 100, 150, 200, 250];
    case 4 % 6 relapses, 28 days each
        days_relapse = 28;
        start_days = [1, 50, 100, 150, 200, 250];
end

timestep_equiv = days_relapse * steps_per_day;
all_windows = [];

for sd = start_days
    start_step = sd * steps_per_day;
    all_windows = [all_windows, start_step:(start_step + timestep_equiv)];
end

%% Determine Peak Probability
if ismember(nn, all_windows)
    peak_prob = 0.0025;  % during relapse period
else
    peak_prob = 0.0002;  % baseline
end

%% Generate New C1 Cells
ypos = 1:299;

rand_nums1 = rand(1, length(ypos));
rand_nums2 = rand(1, length(ypos));

indices1 = rand_nums1 < peak_prob;
indices2 = rand_nums2 < peak_prob;

C1.x = [C1.x, ones(1, sum(indices1)), 2*ones(1, sum(indices2))];
C1.y = [C1.y, ypos(indices1), ypos(indices2)];

num_new_cells = sum(indices1) + sum(indices2);
C1.age = [C1.age, zeros(1, num_new_cells)];

end