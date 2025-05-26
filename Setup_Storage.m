function [spatial] = Setup_Storage(step_cap,C1,C2,C3,myelin)

% Cell 1
spatial.C1x = cell([1 step_cap + 1]); % cell array for C1 x pos 
spatial.C1y = cell([1 step_cap + 1]); % cell array for C1 y pos 
spatial.C1x{1} = C1.x; spatial.C1y{1} = C1.y; % initialise C1 x,y pos 

% Cell 2
spatial.C2x = cell([1 step_cap + 1]); % cell array for C2 x pos 
spatial.C2y = cell([1 step_cap + 1]); % cell array for C2 y pos 
spatial.C2x{1} = C2.x; spatial.C2y{1} = C2.y; % initialise for C2 x,y pos 

% Cell 3 
spatial.C3x = cell([1 step_cap + 1]); % cell array for C3 x pos 
spatial.C3y = cell([1 step_cap + 1]); % cell array for C3 y pos 
spatial.C3x{1} = C3.x; spatial.C3y{1} = C3.y; % initialise for C3 x,y pos

% Myelin 
spatial.myelinstates = cell([1 step_cap + 1]); % cell array for myelin states 
spatial.myelinstates{1} = myelin.state; 

spatial.healthy_x = cell([1 step_cap]); 
spatial.healthy_y = cell([1 step_cap]); 
spatial.intermed_x = cell([1 step_cap]); 
spatial.intermed_y = cell([1 step_cap]); 
spatial.destroy_x = cell([1 step_cap]); 
spatial.destroy_y = cell([1 step_cap]);

%
 spatial.oligo_states = cell([1 step_cap]); 
 spatial.myelin_states = cell([1 step_cap]);


