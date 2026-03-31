% This script runs multiple stochastic simulations (pp) of the agent based model and will store either 
% spatial or population level outputs.
% Please note that results are saved at each timestep

%% Simulation settings

save_spatial_data = 0;     % Choose whether to save spatial outputs or population level results

step_cap = 21600;     % 21600 time steps is 300 days for a single realisation 

num_iterations = 3;       % Number of stochastic runs

relapse_schedule = 4;          % Relapse schedule identifier
oligo_restor     = 1;          % Oligodendrocyte restoration switch 


%% Therapy settings 
% Oligodendrocyte apoptosis and myelination stopping
intervention.oligo_apop    = 14; % damage before oligo undergoes apoptosis
intervention.oligo_stop_my = 10; % damage before oligo stops remyelinating

% Restoration of oligodendrocytes
intervention.restore_oligos   = "on"; % to reintroduce oligodendrocytes
intervention.new_oligo_times  = 5760; % introduction timestep 
intervention.restore_proportion = 1; % proportion of oligodendrocytes to restore (that have undergone apoptosis)

% Oligodendrocyte property changes
intervention.oligo_properties = "on"; % to alter oligodendrocyte properties 
intervention.oligo_prop_times = 5760; % alteration timestep 
intervention.new_apop         = 24; % new apoptosis threshold 
intervention.new_stopmy       = 20; % new stop myelination threshold 

% Blood-brain barrier (BBB) intervention
intervention.BBB           = "on"; % to change BBB permeability 
intervention.BBB_times     = 5760; % alteration timestep 
intervention.new_BBB_prob  = 0.025; % new permeability 

%% Main simulation loop 
for pp = 1:num_iterations
    
    rng(pp);  % Set random seed for reproducibility
    
    % Run simulation
    [result, tracked, myelin, spatial] = main_function( ...
        save_spatial_data, intervention, relapse_schedule, step_cap);
    
    % Store outputs
    bulk_result_data{pp}  = result; % stores cell population data 
    bulk_tracked_data{pp} = tracked; % stores information about tracked oligodendrocytes
    bulk_myelin_data{pp}  = myelin; % stores myelin data 
    
    if save_spatial_data == 1
        bulk_spatial_data{pp} = spatial; % saves spatial coordinates 
    else
        bulk_spatial_data = [];
    end
    
    % Save progress after each iteration (robust to crashes)
    save('test_run.mat', ...
        'bulk_result_data', ...
        'bulk_tracked_data', ...
        'bulk_myelin_data', ...
        'bulk_spatial_data', ...
        'pp', ...
        '-v7.3');
end
