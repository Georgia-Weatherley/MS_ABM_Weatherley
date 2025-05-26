function [result,tracked,myelin,spatial]  = main_script_function(save_spatial,intervention)

threshold_oligo_apop = intervention.oligo_apop; 
threshold_stop_myelination = intervention.oligo_stop_my;
%% Sim Conditions
 step_cap = 21600;     %21600 is 300 days

%%  Setup Simulation
[funcs,PB,sizing,borders,myelin] = Setup_Simulation;
[myelin,C1,C2,C3,tau] = Setup_Agents(myelin,borders,sizing);

%% Setup tracking of individual myelin and oligos 

% num_samples = 2; 
% rowbounds = floor(linspace(min(myelin.yarray),max(myelin.yarray),10));
% colbounds = floor(linspace(min(myelin.xarray),max(myelin.xarray),5));
% 
% rowpos = [];
% colpos = [];
% for ii = 1:length(rowbounds)-1
%     for jj = 1:length(colbounds)-1
%         rowpos = [rowpos randi([rowbounds(ii) rowbounds(ii+1)],1,num_samples)];
%         colpos = [colpos randi([colbounds(jj) colbounds(jj+1)],1,num_samples)];
%     end
% end
% 
% tracked.indices = [];
% for ff = 1:length(rowpos)
%     tracked.indices = [tracked.indices find(myelin.xarray == colpos(ff) & myelin.yarray == rowpos(ff))];
%     tracked.oligo_indices = myelin.oligo_tracker(tracked.indices);
% end

indx1 = [977 935 1036];
indx2 = [133 168 232]; 
indx3 = [429 486 529];

m_indx1 = [24413 23363 25888];
m_indx2 = [3313 4188 5788];
m_indx3 = [10713 12138 13213];

tracked.indices = [m_indx1 m_indx2 m_indx3];
tracked.oligo_indices = [indx1 indx2 indx3];


%% Setup Storage: 
tracked.states = zeros(length(tracked.indices),step_cap); 
tracked.oligo_states = zeros(length(tracked.oligo_indices),step_cap); 

result.myelinating_oligo = zeros(1,step_cap); 
result.dead_oligo = zeros(1,step_cap); 
result.nonmyelinating_oligo = zeros(1,step_cap); 
result.healthy_myelin = zeros(1,step_cap); 
result.intermediate_myelin = zeros(1,step_cap); 
result.destroyed_myelin = zeros(1,step_cap); 

result.C1_pop = zeros(1,step_cap); 
result.C2_pop = zeros(1,step_cap); 
result.C3_pop = zeros(1,step_cap); 

if save_spatial == true 
    [spatial] = Setup_Storage(step_cap,C1,C2,C3,myelin);
else 
    spatial = 0;
end
%% Run each iteration  
tic
 for nn = 1:step_cap 
     nn
     if intervention.active == "on"
         if intervention.BBB == "on" && (sum(nn == intervention.BBB_times) == 1)
             PB.leavingprob = intervention.new_BBB_prob; 
         end 

         if intervention.restore_oligos == "on" && (sum(nn == intervention.new_oligo_times) == 1)
             myelin = new_oligos(myelin);
         end

         if intervention.oligo_properties == "on" && (sum(nn == intervention.oligo_prop_times) == 1)
            threshold_oligo_apop = intervention.new_apop;
            threshold_stop_myelination = intervention.new_stopmy;    
         end
     end
    
    % Move C3 with unbiased or biased movement 
    if length(C3.x) ~= 0
        myelin.sum = sum(myelin.state);
        if myelin.sum == myelin.fullyhealthyquant || myelin.sum == 0
            [C3.x,C3.y]= UNBIASED(C3.x,C3.y,PB);
        else
            C3 = BIASED(C3,myelin,funcs,PB);
        end
        
    end
    % Move C1 and C2 cells in unbiased manner 
    [C1.x,C1.y] = UNBIASED(C1.x,C1.y,PB); 
    [C2.x,C2.y]= UNBIASED(C2.x,C2.y,PB); 
    wrap_down.C1 = C1.y == sizing.domainheight; C1.y(wrap_down.C1) = 1; 
    wrap_down.C2 = C2.y == sizing.domainheight; C2.y(wrap_down.C2) = 1; 
    wrap_down.C3 = C3.y == sizing.domainheight; C3.y(wrap_down.C3) = 1; 
    wrap_up.C1 = C1.y == 0; C1.y(wrap_up.C1) = sizing.domainheight-1; 
    wrap_up.C2 = C2.y == 0; C2.y(wrap_up.C2) = sizing.domainheight-1; 
    wrap_up.C3 = C3.y == 0; C3.y(wrap_up.C3) = sizing.domainheight-1; 

    C1.x(C1.x == 0) = 1; % reflecting boundary on LHS 

    C2.x(C2.x == borders.rightcns+1) = borders.rightcns; % reflecting boundary on RHS

    evict.C1 = C1.x == sizing.domainwidth| C1.x == 0; % mask to enforce the boundaries
    evict.C3 = C3.x == sizing.domainwidth;  
   
    C1.x = C1.x(~evict.C1);  C1.y = C1.y(~evict.C1); % Remove exited cells 
    C3.x = C3.x(~evict.C3);  C3.y = C3.y(~evict.C3);

    myelin = REMYELINATION_WITH_OLIGO(threshold_stop_myelination,myelin,funcs) ;
    myelin = OLIGO_APOPTOSIS(myelin, threshold_oligo_apop,threshold_stop_myelination);
    
    if length(C3.x) ~= 0
        myelin = DEMYELINATION(myelin,C3); % Check for myelin degeneration 
    end
    
    if length(C1.x) ~= 0 & length(C2.x) ~= 0 % potential to create C3 cells 
        C1.pos = [C1.x;C1.y]'; % Col 1 is x pos, Col 2 is y pos of C1
        C2.pos = [C2.x;C2.y]'; % Col 1 is x pos, Col 2 is y pos of C2 
        [~,index_C1,~] = intersect(C1.pos,C2.pos,'rows'); % search rows to look for identical x,y pos of C1 and C2   
        new_C3x_pos = C1.x(index_C1); % create a C3 cell at the intersection point   
        C3.x = [C3.x new_C3x_pos]; % add x coords to C3 
        C3.y = [C3.y C1.y(index_C1)]; % add y coords to C3  
        C1.x = C1.x(setdiff(1:end,index_C1)); % remove x coords from C1
        C1.y = C1.y(setdiff(1:end,index_C1)); % remove y coords from C1
    end

    C1 = SITE_METHOD(C1,nn,tau);

    mort.C1mask = rand(1,length(C1.x)) < C1.deathrate; % Create logical masks of dead cells 
    mort.C3mask = rand(1,length(C3.x)) < C3.deathrate; 
   
    C1.x = C1.x(~mort.C1mask); C1.y = C1.y(~mort.C1mask);  
    C3.x = C3.x(~mort.C3mask); C3.y = C3.y(~mort.C3mask); 
    
    % save results
    result.myelinating_oligo(nn) = sum(myelin.oligo_state == 1);
    result.dead_oligo(nn) = sum(myelin.oligo_state == 0);
    result.nonmyelinating_oligo(nn) = sum(myelin.oligo_state == 2);

    result.healthy_myelin(nn) = sum(myelin.state == myelin.grades); 
    result.intermediate_myelin(nn) = sum(myelin.state > 0 & myelin.state < myelin.grades);
    result.destroyed_myelin(nn) = sum(myelin.state == 0);

    result.C1_pop(nn) = length(C1.x); 
    result.C2_pop(nn) = length(C2.x); 
    result.C3_pop(nn) = length(C3.x); 

    tracked.states(:,nn) = myelin.state(tracked.indices);
    tracked.oligo_states(:,nn) = myelin.oligo_state(tracked.oligo_indices);

    if save_spatial == true 

    spatial.C1x{nn+1} = C1.x; spatial.C1y{nn+1} = C1.y; % cell array for C1 x,y pos 
    spatial.C2x{nn+1} = C2.x; spatial.C2y{nn+1} = C2.y; % cell array for C2 x,y pos 
    spatial.C3x{nn+1} = C3.x; spatial.C3y{nn+1} = C3.y; % cell array for C3 x,y pos

    spatial.healthy_x{nn} = myelin.xarray(myelin.state == myelin.grades);
    spatial.healthy_y{nn} = myelin.yarray(myelin.state == myelin.grades);
    spatial.intermed_x{nn} = myelin.xarray(myelin.state < myelin.grades & myelin.state > 0);
    spatial.intermed_y{nn} = myelin.yarray(myelin.state < myelin.grades & myelin.state > 0);
    spatial.destroy_x{nn} = myelin.xarray(myelin.state == 0);
    spatial.destroy_y{nn} = myelin.yarray(myelin.state== 0);

    spatial.oligo_states{nn} = myelin.oligo_state;
    spatial.myelin_states{nn} = myelin.state;
    end
 end

 if save_spatial == true 
     [spatial] = tidy_to_plot(spatial); 
 end
