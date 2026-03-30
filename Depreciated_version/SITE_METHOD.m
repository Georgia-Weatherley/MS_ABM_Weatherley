function [C1] = SITE_METHOD(C1,nn,tau)
%SITE_METHOD Introduces new C1 cells (T cells) in the peripheral blood. 

mins_per_day = 1440;
steps_per_day = mins_per_day/tau; % number of timesteps in one day
 
num_weeks = 4; % weeks of relapse 
timestep_equiv = num_weeks*7*mins_per_day /tau; 

start1 = 1; % day when relapse begins 
start2 = 100; % day when relapse begins 
start3 = 200; % day when relapse begins 

window1 = start1:1:timestep_equiv;
window2 =  start2*steps_per_day:1: (mins_per_day*start2/ tau)+timestep_equiv; 
window3 = start3*steps_per_day:1:(mins_per_day*start3/ tau) +timestep_equiv;

all_windows = [window1 window2 window3]; 

if ismember(nn,all_windows) % it is a relapsing period 
%      peak_prob = 0.0027;
     peak_prob = 0.0025;
     rand_nums1 = rand(1,299);
     rand_nums2 = rand(1,299);
     indices1 = rand_nums1 < peak_prob; 
     indices2 = rand_nums2 < peak_prob; 
     ypos = 1:1:299; 
     C1.x = [C1.x ones(1,sum(indices1)) 2*ones(1,sum(indices2))];
     C1.y = [C1.y ypos(indices1) ypos(indices2)]; 
     num_new_cells = sum(indices1) + sum(indices2);
     C1.age = [C1.age zeros(1,num_new_cells)]; 
else % it is not a relapsing period 
     peak_prob = 0.0002; 
     rand_nums1 = rand(1,299);
     rand_nums2 = rand(1,299);
     indices1 = rand_nums1 < peak_prob; 
     indices2 = rand_nums2 < peak_prob; 
     ypos = 1:1:299; 
     C1.x = [C1.x ones(1,sum(indices1)) 2*ones(1,sum(indices2))];
     C1.y = [C1.y ypos(indices1) ypos(indices2)]; 
     num_new_cells = sum(indices1) + sum(indices2);
     C1.age = [C1.age zeros(1,num_new_cells)];  
end


