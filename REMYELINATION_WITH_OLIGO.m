function [myelin] = REMYELINATION_WITH_OLIGO(local_threshold,myelin,funcs) 

reshaped_myelin = reshape(myelin.state, myelin.oligo_dim^2, myelin.oligo_counter);
indicate_destroyed_myelin = reshaped_myelin == 0;
total_damage_to_oligo = sum(indicate_destroyed_myelin); 

% matrix where cols hold damage specific to each oligo
dim_fix_damage = repmat(total_damage_to_oligo,myelin.oligo_dim^2,1); 

% find which oligos are healthy
oligo_healthy = myelin.oligo_state == 1; 
dim_fix_oligo = repmat(oligo_healthy,myelin.oligo_dim^2,1); 
array_healthy_oligos = reshape(dim_fix_oligo,1,[]);

% if enough neighbours, increase time waiting for regen (increase prob of
% regen) 
myelin.timer(array_healthy_oligos) = myelin.timer(array_healthy_oligos) + 1;

% Generate random numbers to sample regeneration function
% timer_rand = rand(1,length(myelin.state)); 
% regen_func_eval = funcs.myelin(myelin.timer); % values of regen func evaluated at timers 
% timer_mask = timer_rand < regen_func_eval; % see if under func threshold 

timer_mask = myelin.timer >= myelin.healtime;

% matrix where rows hold time specific to each myelin, cols are sorted by
% oligo
dim_fix_timer = reshape(timer_mask, myelin.oligo_dim^2, myelin.oligo_counter);

build_back_mask = dim_fix_damage < local_threshold & dim_fix_timer == 1 & dim_fix_oligo == 1;
array_build_back_mask = reshape(build_back_mask,1,[]); 

% restore the healthy myelin
myelin.state(array_build_back_mask) = myelin.state(array_build_back_mask) + 1;

% Reset the regeneration timer of myelin just healed 
myelin.timer(array_build_back_mask) = 0; 

% Stop myelin from over healing
myelin.state(myelin.state > myelin.grades) = myelin.grades;

% Take myelin that is now fully healthy out of the regeneration timer  
myelin.timer(myelin.state == myelin.grades) = -1; 