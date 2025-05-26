
%% #1: Run for: general, untreated case 
save_spatial_data = false; 

intervention.active = "on";

intervention.oligo_apop = 14;
intervention.oligo_stop_my = 10; 

intervention.restore_oligos = "off";
intervention.new_oligo_times = 5760;

intervention.oligo_properties = "off";
intervention.oligo_prop_times = 5760;
intervention.new_apop = 24;
intervention.new_stopmy = 20;

intervention.BBB = "off"; 
intervention.BBB_times = 5760;
intervention.new_BBB_prob = 0.025;

num_iterations = 1; 

for pp = 1:num_iterations
        rng(10)
        if save_spatial_data == false
        [result,tracked,myelin,spatial]  = main_script_function(save_spatial_data,intervention);
            bulk_result_data{pp} = result;
            bulk_tracked_data{pp} = tracked;
            bulk_myelin_data{pp} = myelin;
            bulk_spatial_data = []; 
        else
            [result,tracked,myelin,spatial]  = main_script_function(save_spatial_data,intervention);
            bulk_result_data{pp} = result;
            bulk_tracked_data{pp} = tracked;
            bulk_myelin_data{pp} = myelin;
            bulk_spatial_data{pp} = spatial; 
        end
       
     save('one_iteration_14_10.mat','bulk_result_data','bulk_tracked_data','bulk_myelin_data','bulk_spatial_data','pp','-v7.3')
     pp
end
   