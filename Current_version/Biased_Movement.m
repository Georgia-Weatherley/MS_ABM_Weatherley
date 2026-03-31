function [C3] = Biased_Movement(C3, myelin, funcs, PB)

% BIASED_MOVEMENT
%
% Performs biased movement of C3 agents toward nearby healthy myelin.
% Bias strength depends on distance and directional components.
%
% INPUTS:
%   C3     - structure containing C3 agent positions
%   myelin - myelin structure (state and spatial information)
%   funcs  - structure containing bias functions
%   PB     - probabilistic boundary parameters
%
% OUTPUT:
%   C3 - updated C3 structure with new positions

%% Extract Data
C3_x = C3.x;
C3_y = C3.y;

C3_bias_func = funcs.C3_bias;

delta = 1;   % Step size

% Combine coordinates
C3_xy = [C3_x; C3_y]';

%% Identify Myelin 
alive_my_x = myelin.xarray(myelin.state > 0);
alive_my_y = myelin.yarray(myelin.state > 0);

living_my_xy = [alive_my_x; alive_my_y]';

[length_C3_array, ~] = size(C3_xy);

%% Distance to Nearest Myelin
dist    = zeros(1, length_C3_array);
k_index = zeros(1, length_C3_array);

for rr = 1:length_C3_array
    
    cell_of_interest = C3_xy(rr, :);
    
    fill      = (living_my_xy - cell_of_interest)';
    norm_dist = sqrt(sum(fill.^2));
    
    % Identify nearest myelin (random tie-breaking)
    myel_index_options = find(norm_dist == min(norm_dist));
    choice = randi(length(myel_index_options), 1);
    
    k_index(rr) = myel_index_options(choice);
    dist(rr)    = norm_dist(k_index(rr));
end

%% Directional Differences
x_differences = alive_my_x(k_index) - C3_x;
y_differences = alive_my_y(k_index) - C3_y;

% Direction masks
bias_east  = x_differences >= 0;
bias_west  = x_differences < 0;
bias_north = y_differences >= 0;
bias_south = y_differences < 0;

%% Bias Weighting
x_bias_ratio = abs(x_differences) ./ ...
               (abs(x_differences) + abs(y_differences));

y_bias_ratio = abs(y_differences) ./ ...
               (abs(x_differences) + abs(y_differences));

% Distance-based bias
prob_to_distribute = C3_bias_func(dist);

% Base probability
initial_bias = 0.25 * ones(1, length(C3_x));

% Additional bias components
x_extra_bias = x_bias_ratio .* prob_to_distribute;
y_extra_bias = y_bias_ratio .* prob_to_distribute;

% Total directional bias
sum_x_bias = initial_bias + x_extra_bias;
sum_y_bias = initial_bias + y_extra_bias;

%% Distance Threshold
max_dist = 80;

close_enough = dist' <= max_dist;

%% Assign Directional Probabilities
bias_east_amount  = bias_east  .* sum_x_bias .* close_enough';
bias_west_amount  = bias_west  .* sum_x_bias .* close_enough';
bias_north_amount = bias_north .* sum_y_bias .* close_enough';
bias_south_amount = bias_south .* sum_y_bias .* close_enough';

%% Redistribute Remaining Probability
quan_no_bias = (bias_east_amount  == 0) + ...
               (bias_west_amount  == 0) + ...
               (bias_north_amount == 0) + ...
               (bias_south_amount == 0);

fill_in_right_mask = bias_east_amount  == 0;
fill_in_left_mask  = bias_west_amount  == 0;
fill_in_up_mask    = bias_north_amount == 0;
fill_in_down_mask  = bias_south_amount == 0;

bias_distributed = bias_east_amount + bias_west_amount + ...
                   bias_north_amount + bias_south_amount;

leftover_probs = 1 - bias_distributed;

% Evenly distribute remaining probability
bias_east_amount(fill_in_right_mask) = ...
    leftover_probs(fill_in_right_mask) ./ ...
    quan_no_bias(fill_in_right_mask);

bias_west_amount(fill_in_left_mask) = ...
    leftover_probs(fill_in_left_mask) ./ ...
    quan_no_bias(fill_in_left_mask);

bias_north_amount(fill_in_up_mask) = ...
    leftover_probs(fill_in_up_mask) ./ ...
    quan_no_bias(fill_in_up_mask);

bias_south_amount(fill_in_down_mask) = ...
    leftover_probs(fill_in_down_mask) ./ ...
    quan_no_bias(fill_in_down_mask);

%% Cumulative Probabilities
W_n_E_sum     = bias_west_amount + bias_east_amount;
W_n_E_n_S_sum = W_n_E_sum + bias_south_amount;

%% Sample Movement
r_C3 = rand(1, length(C3_x));

% Direction masks
C3_west_mask  = r_C3 < bias_west_amount;

if sum(C3_west_mask) ~= 0
    r_BBB_west = rand(1, length_C3_array);
    accept_mask = r_BBB_west <= PB.enteringprob;
    
    C3_west_mask(C3_west_mask == 1 & ...
                 C3_x == PB.rightside & ...
                 accept_mask ~= 1) = 0;
end

C3_east_mask = r_C3 >= bias_west_amount & r_C3 < W_n_E_sum;

if sum(C3_east_mask) ~= 0
    r_BBB_east = rand(1, length_C3_array);
    accept_mask = r_BBB_east <= PB.leavingprob;
    
    C3_east_mask(C3_east_mask == 1 & ...
                 C3_x == PB.leftside & ...
                 accept_mask ~= 1) = 0;
end

C3_south_mask = r_C3 >= W_n_E_sum & r_C3 < W_n_E_n_S_sum;
C3_north_mask = r_C3 >= W_n_E_n_S_sum;

%% Apply Movement
C3_x(C3_west_mask)  = C3_x(C3_west_mask)  - delta;
C3_x(C3_east_mask)  = C3_x(C3_east_mask)  + delta;
C3_y(C3_south_mask) = C3_y(C3_south_mask) - delta;
C3_y(C3_north_mask) = C3_y(C3_north_mask) + delta;

%% Update Output
C3.x = C3_x;
C3.y = C3_y;

end