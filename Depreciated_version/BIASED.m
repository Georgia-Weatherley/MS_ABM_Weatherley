function [C3] = BIASED(C3,myelin,funcs,PB)
C3_x = C3.x; 
C3_y = C3.y;
C3_bias_func = funcs.C3_bias;
% define step size 
delta = 1; 

% Create 2 column matrix where column 1 is x pos, column 2 is y pos   
C3_xy = [C3_x;C3_y]'; 

% Create a mask of x,y pos where myelin is not dead 
alive_my_x = myelin.xarray(myelin.state> 0); % array created using mask 
alive_my_y = myelin.yarray(myelin.state> 0); % array created using mask 

% Create 2 column matrix where column 1 is x pos, column 2 is y pos 
living_my_xy = [alive_my_x;alive_my_y]';

[length_C3_array,~] = size(C3_xy);

dist = zeros(1,length_C3_array); 
k_index = zeros(1,length_C3_array); 

for rr = 1:length_C3_array
    cell_of_interest = C3_xy(rr,:);  
    fill = (living_my_xy - cell_of_interest)';
    norm_dist = sqrt(sum(fill.^2));
%     norm_dist = vecnorm(fill);  % replace with formula
  %  [~,myelin_index_options] = min(norm_dist);
  myel_index_options = find(norm_dist == min(norm_dist));
    choice = randi(length(myel_index_options),1)  ;
    
    k_index(rr) = myel_index_options(choice);    
    dist(rr) = norm_dist(k_index(rr)) ;
end 

% Take the x and y distance between the C3 cells and their closest myelin
x_differences = alive_my_x(k_index) - C3_x; % take difference of cell and myelin x pos
y_differences = alive_my_y(k_index) - C3_y; % take difference of cell and myelin y pos 


% Create logical masks of bias directions 
bias_east = x_differences >= 0; 
bias_west = x_differences < 0; 
bias_north = y_differences >= 0; 
bias_south = y_differences < 0;

% Work out how much x,y components contribute to distance
x_bias_ratio = abs(x_differences)./(abs(x_differences)+abs(y_differences)); 
y_bias_ratio = abs(y_differences)./(abs(x_differences)+abs(y_differences));

% Sample bias function to determine amount of bias based on distance 
prob_to_distribute = C3_bias_func(dist); 

% Create baseline prob of 25% chance in each direction 
initial_bias = 0.25*ones(1,length(C3_x));

% Work out amount of additional bias based on x,y component ratios and
% distance based bias amount 
x_extra_bias = x_bias_ratio.*prob_to_distribute ;
y_extra_bias = y_bias_ratio.*prob_to_distribute ;

% total bias in each direction 
sum_x_bias = initial_bias + x_extra_bias ; 
sum_y_bias = initial_bias + y_extra_bias ; 

% Set the max distance from which myelin can bias C3 cells   
max_dist = 80;

% Create logical mask of C3 cells close enough to myelin to be biased but
% remove cells that I want to keep in the same spot
close_enough = dist' <= max_dist;

% split bias into each direction and only apply if C3 cells close enough
% 1st term logical, 2nd term a quantity, 3rd term a logical 
bias_east_amount = bias_east.*sum_x_bias.*close_enough' ; 
bias_west_amount = bias_west.*sum_x_bias.*close_enough' ; 
bias_north_amount = bias_north.*sum_y_bias.*close_enough' ; 
bias_south_amount = bias_south.*sum_y_bias.*close_enough' ; 

% Find how many directions unaffected by myelin bias 
quan_no_bias = (bias_east_amount == 0)+(bias_west_amount == 0) + ...
                   (bias_north_amount == 0) + (bias_south_amount == 0); 
               
% Find where these unaffected/unbiased directions are 
fill_in_right_mask = bias_east_amount == 0; 
fill_in_left_mask = bias_west_amount == 0; 
fill_in_up_mask = bias_north_amount == 0; 
fill_in_down_mask = bias_south_amount == 0;   

% Find amount of bias already distributed 
bias_distributed = bias_east_amount + bias_west_amount + ... 
                        bias_north_amount + bias_south_amount; 
                    
% Work out remaining RW prob left to distribute to other directions 
leftover_probs = 1 - bias_distributed ;      

 % Fill in unassigned directions with equal split of leftover prob 
bias_east_amount(fill_in_right_mask) =  leftover_probs(fill_in_right_mask)./quan_no_bias(fill_in_right_mask); 
bias_west_amount(fill_in_left_mask) =  leftover_probs(fill_in_left_mask)./quan_no_bias(fill_in_left_mask); 
bias_north_amount(fill_in_up_mask) =  leftover_probs(fill_in_up_mask)./quan_no_bias(fill_in_up_mask); 
bias_south_amount(fill_in_down_mask) =  leftover_probs(fill_in_down_mask)./quan_no_bias(fill_in_down_mask);            
 
% Create brackets 
W_n_E_sum = bias_west_amount + bias_east_amount; 
W_n_E_n_S_sum = W_n_E_sum + bias_south_amount; 

% Generate random probabilities   
r_C3 = rand(1,length(C3_x)); % horizontal array

% Create mask for each direction 
C3_west_mask = r_C3 < bias_west_amount; % determine which to shift left  
if sum(C3_west_mask) ~= 0 
    r_BBB_west = rand(1,length_C3_array);
    accept_mask = r_BBB_west <= PB.enteringprob;
    C3_west_mask(C3_west_mask == 1 & C3_x == PB.rightside & accept_mask ~= 1) = 0;
end

C3_east_mask = r_C3 >= bias_west_amount & r_C3 < W_n_E_sum; % determine which to shift right  
if sum(C3_east_mask) ~= 0 
    r_BBB_east = rand(1,length_C3_array);
    accept_mask = r_BBB_east <= PB.leavingprob;
    C3_east_mask(C3_east_mask == 1 & C3_x == PB.leftside & accept_mask ~= 1) = 0;
end

C3_south_mask = r_C3 >=  W_n_E_sum & r_C3 < W_n_E_n_S_sum; % determine which to shift down            
C3_north_mask = r_C3 >= W_n_E_n_S_sum; % determine which to shift up  

% Shift C3 cells 
C3_x(C3_west_mask)  = C3_x(C3_west_mask)  - delta; % shift left
C3_x(C3_east_mask)  = C3_x(C3_east_mask)  + delta; % shift right
C3_y(C3_south_mask) = C3_y(C3_south_mask) - delta; % shift down 
C3_y(C3_north_mask) = C3_y(C3_north_mask) + delta; % shift up 

C3.x = C3_x; 
C3.y = C3_y;


end 