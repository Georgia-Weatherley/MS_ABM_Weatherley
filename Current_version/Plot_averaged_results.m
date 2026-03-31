% This script is used to plot population-level results of the ABM averaged
% across multiple realisations. 
% It requites data from the ABM_main_script to be in the workspace, or
% loaded in as a .mat file

%% Confirm Simulation Parameters 
num_iter = length(bulk_result_data); % number of realisations averaged over 
step_cap = length(bulk_result_data{1}.C1_pop); % number of time steps used in simulation 
tau = 20; % time step duration 
mins_per_day = 1440;
steps_per_day = mins_per_day/tau; % number of timesteps in one day

relapse_schedule = 4; % relapse schedule used (timing and duration)

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

timestep_equiv = days_relapse*mins_per_day /tau; 

%% Set up the background patch to indicate relapse on graph 
num_relapses = length(start_days);
Patch_X = zeros(num_relapses, 4); % preallocate

for i = 1:num_relapses
    start_step = start_days(i) * steps_per_day; % start in timesteps
    OG_window = start_step : (start_step + timestep_equiv);
    
    % Create the patch X coordinates for this relapse
    Patch_X(i,:) = [OG_window(1), OG_window(end), OG_window(end), OG_window(1)] / steps_per_day;
end

% Y coordinates are the same for all patches
Patch_Y = [-1 -1 1000 1000];

%% Set up empty matrices to collate the data
cell_1_mat = zeros(num_iter,step_cap);
cell_2_mat =zeros(num_iter,step_cap);
cell_3_mat = zeros(num_iter,step_cap);
myel_healthy_mat = zeros(num_iter,step_cap);
myel_intermed_mat = zeros(num_iter,step_cap);
myel_destroy_mat = zeros(num_iter,step_cap);
oligo_dead_mat = zeros(num_iter,step_cap);
oligo_myel_mat = zeros(num_iter,step_cap);
oligo_nonmyel_mat = zeros(num_iter,step_cap);

min_C3 = 1000*ones(1,step_cap); 
max_C3 = zeros(1,step_cap); 
min_C1 = 1000*ones(1,step_cap); 
max_C1 = zeros(1,step_cap); 
min_C2 = 1000*ones(1,step_cap); 
max_C2 = zeros(1,step_cap); 
min_healthy_myel = 10000000*ones(1,step_cap); 
max_healthy_myel = zeros(1,step_cap); 
min_itermed_myel = 100000000*ones(1,step_cap); 
max_itermed_myel = zeros(1,step_cap); 
min_destroy_myel = 100000000*ones(1,step_cap); 
max_destroy_myel = zeros(1,step_cap); 

min_oligo_dead = 10000000*ones(1,step_cap); 
max_oligo_dead = zeros(1,step_cap); 
min_oligo_myel = 10000000*ones(1,step_cap); 
max_oligo_myel = zeros(1,step_cap); 
min_oligo_nonmyel = 100000000*ones(1,step_cap); 
max_oligo_nonmyel = zeros(1,step_cap); 

%% Extract the simulation data to be compatible with plotting 
for ii = 1:num_iter 
    cell_3_mat(ii,:) = bulk_result_data{1,ii}.C3_pop;
    min_C3 = min(min_C3,bulk_result_data{1,ii}.C3_pop);
    max_C3 = max(max_C3,bulk_result_data{1,ii}.C3_pop);

    cell_1_mat(ii,:) = bulk_result_data{1,ii}.C1_pop;
    min_C1 = min(min_C1,bulk_result_data{1,ii}.C1_pop);
    max_C1 = max(max_C1,bulk_result_data{1,ii}.C1_pop);

    cell_2_mat(ii,:) = bulk_result_data{1,ii}.C2_pop;
    min_C2 = min(min_C2,bulk_result_data{1,ii}.C2_pop);
    max_C2 = max(max_C2,bulk_result_data{1,ii}.C2_pop);

    myel_healthy_mat(ii,:) = bulk_result_data{1,ii}.healthy_myelin;
    min_healthy_myel = min(min_healthy_myel,bulk_result_data{1,ii}.healthy_myelin);
    max_healthy_myel = max(max_healthy_myel,bulk_result_data{1,ii}.healthy_myelin);

    myel_intermed_mat(ii,:) = bulk_result_data{1,ii}.intermediate_myelin;
    min_itermed_myel = min(min_itermed_myel,bulk_result_data{1,ii}.intermediate_myelin);
    max_itermed_myel = max(max_itermed_myel,bulk_result_data{1,ii}.intermediate_myelin);

    myel_destroy_mat(ii,:) = bulk_result_data{1,ii}.destroyed_myelin;
    min_destroy_myel = min(min_destroy_myel,bulk_result_data{1,ii}.destroyed_myelin);
    max_destroy_myel = max(max_destroy_myel,bulk_result_data{1,ii}.destroyed_myelin);

    oligo_dead_mat(ii,:) = bulk_result_data{1,ii}.dead_oligo;
    min_oligo_dead = min(min_oligo_dead,bulk_result_data{1,ii}.dead_oligo);
    max_oligo_dead = max(max_oligo_dead,bulk_result_data{1,ii}.dead_oligo);

    oligo_myel_mat(ii,:) = bulk_result_data{1,ii}.myelinating_oligo; 
    min_oligo_myel = min(min_oligo_myel,bulk_result_data{1,ii}.myelinating_oligo);
    max_oligo_myel = max(max_oligo_myel,bulk_result_data{1,ii}.myelinating_oligo);


    oligo_nonmyel_mat(ii,:) = bulk_result_data{1,ii}.nonmyelinating_oligo; 
    min_oligo_nonmyel = min(min_oligo_nonmyel,bulk_result_data{1,ii}.nonmyelinating_oligo);
    max_oligo_nonmyel = max(max_oligo_nonmyel,bulk_result_data{1,ii}.nonmyelinating_oligo);
end

%% Average the data 
divisor_array = (1:num_iter);

cell_3_av = cumsum(cell_3_mat);
cell_3_av = cell_3_av./divisor_array';

cell_1_av = cumsum(cell_1_mat);
cell_1_av = cell_1_av./divisor_array';

cell_2_av = cumsum(cell_2_mat);
cell_2_av = cell_2_av./divisor_array';

myel_healthy_av = cumsum(myel_healthy_mat); 
myel_healthy_av  = myel_healthy_av./divisor_array';

myel_intermed_av = cumsum(myel_intermed_mat); 
myel_intermed_av  = myel_intermed_av./divisor_array';

myel_destroy_av = cumsum(myel_destroy_mat); 
myel_destroy_av  = myel_destroy_av./divisor_array';

oligo_dead_av = cumsum(oligo_dead_mat); 
oligo_dead_av  = oligo_dead_av./divisor_array';

oligo_nonmyel_av = cumsum(oligo_nonmyel_mat); 
oligo_nonmyel_av  = oligo_nonmyel_av./divisor_array';

oligo_myel_av = cumsum(oligo_myel_mat); 
oligo_myel_av  = oligo_myel_av./divisor_array';

%% Plot Averaged Cell Populations 
[funcs,PB,sizing,borders,myelin] = Setup_Simulation;
[col,patching,plotstyle] = colours_and_patches(borders,sizing);

figure, hold on 
% Shade relapse windows
for ii = 1:num_relapses
    p1 = patch(Patch_X(ii,:),Patch_Y,'k','FaceAlpha',0.09,'LineStyle','none');
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
x = (1:step_cap)/(1440/20);  
% Plot C1 data 
y =  cell_1_av(end,:);
plot(x,y','Color',col.C1,'DisplayName','T cell')  
f = fill([x fliplr(x)],[max_C1 fliplr(min_C1)],1,'FaceColor', col.C1, 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot C2 data 
y =  cell_2_av(end,:);
plot(x,y','Color',col.C2,'DisplayName','Macrophage')
f = fill([x fliplr(x)],[max_C2 fliplr(min_C2)],1,'FaceColor', col.C2, 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot C3 data 
y =  cell_3_av(end,:);
plot(x,y','k-','Color',col.C3,'DisplayName','Activated T cell')
hold on
f = fill([x fliplr(x)],[max_C3 fliplr(min_C3)],1,'FaceColor', col.C3, 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([0 270]);
xlim([1 step_cap/(1440/20)])
ylabel('Cell count')
xlabel('Time [days]')
legend boxoff   
set(gca,'FontSize',20,"XMinorTick", "on", "YMinorTick", "on", "XMinorGrid","on", "YMinorGrid", "on");

%% Plot Averaged Myelin Populations 
myel_cols = 1/255 * [120,198,121; 65,171,93; 35,132,67; 0,104,55; 0,69,41];
figure, hold on 
% Shade relapse windows
for ii = 1:num_relapses
    p1 = patch(Patch_X(ii,:),Patch_Y,'k','FaceAlpha',0.09,'LineStyle','none');
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
x = (1:step_cap)/(1440/20); 
% Plot fully intact myelin data
y = 100*myel_healthy_av(end,:)/bulk_myelin_data{1}.numpieces;
plot(x,y,'DisplayName','Fully intact','Color',myel_cols(1,:)) 
f = fill([x fliplr(x)],[100*max_healthy_myel/bulk_myelin_data{1}.numpieces 100*fliplr(min_healthy_myel)/bulk_myelin_data{1}.numpieces],1,'FaceColor',myel_cols(1,:), 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot partially intact myelin data
y =  100*myel_intermed_av(end,:)/bulk_myelin_data{1}.numpieces;
plot(x,y,'DisplayName','Partially damaged','Color',myel_cols(3,:))
f = fill([x fliplr(x)],[100*max_itermed_myel/bulk_myelin_data{1}.numpieces 100*fliplr(min_itermed_myel)/bulk_myelin_data{1}.numpieces],1,'FaceColor',myel_cols(3,:), 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot fully degraded myelin data
y = 100*myel_destroy_av(end,:)/bulk_myelin_data{1}.numpieces;
plot(x,y,'DisplayName','Fully damaged','Color',myel_cols(5,:))
f = fill([x fliplr(x)],[100*max_destroy_myel/bulk_myelin_data{1}.numpieces 100*fliplr(min_destroy_myel)/bulk_myelin_data{1}.numpieces],1,'FaceColor',myel_cols(5,:), 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlim([1 step_cap/(1440/20)])
ylabel('Proportion of population')
xlabel('Time [days]')
legend boxoff   
ylim([0 100])
set(gca,'FontSize',20,"XMinorTick", "on", "YMinorTick", "on", "XMinorGrid","on", "YMinorGrid", "on");

%% Plot Averaged Oligo Populations  

figure, hold on 
oligo_cols = 1/255 * [247,104,161; 221,52,151; 174,1,126; 122,1,119; 73,0,106];
% Shade relapse windows
for ii = 1:num_relapses
    p1 = patch(Patch_X(ii,:),Patch_Y,'k','FaceAlpha',0.09,'LineStyle','none');
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
x = (1:step_cap)/(1440/20); 
% Plot myelinating oligo data
y = 100*oligo_myel_av(end,:)/bulk_myelin_data{1}.oligo_counter;
plot(x,y,'DisplayName','Myelinating','Color',oligo_cols(1,:)); 
f = fill([x fliplr(x)],[100*max_oligo_myel/bulk_myelin_data{1}.oligo_counter 100*fliplr(min_oligo_myel)/bulk_myelin_data{1}.oligo_counter],1,'FaceColor',oligo_cols(1,:), 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot nonmyelinating data
y = 100*oligo_nonmyel_av(end,:)/bulk_myelin_data{1}.oligo_counter;
plot(x,y,'DisplayName','Nonmyelinating','Color',oligo_cols(3,:));
f = fill([x fliplr(x)],[100*max_oligo_nonmyel/bulk_myelin_data{1}.oligo_counter 100*fliplr(min_oligo_nonmyel)/bulk_myelin_data{1}.oligo_counter],1,'FaceColor',oligo_cols(3,:), 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot dead oligo data 
y = 100*oligo_dead_av(end,:)/bulk_myelin_data{1}.oligo_counter;
plot(x,y,'DisplayName','Undergone apoptosis','Color',oligo_cols(5,:));
f = fill([x fliplr(x)],[100*max_oligo_dead/bulk_myelin_data{1}.oligo_counter 100*fliplr(min_oligo_dead)/bulk_myelin_data{1}.oligo_counter],1,'FaceColor',oligo_cols(5,:), 'Edgecolor', 'none', 'facealpha', 0.4,'LineStyle','none');
f.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([1 step_cap/(1440/20)])
ylabel('Proportion of population')
legend boxoff  
ylim([0 100])
xlabel('Time [days]')
set(gca,'FontSize',20,"XMinorTick", "on", "YMinorTick", "on", "XMinorGrid","on", "YMinorGrid", "on");
