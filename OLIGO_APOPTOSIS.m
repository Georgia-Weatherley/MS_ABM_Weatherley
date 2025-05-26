function myelin = OLIGO_APOPTOSIS(myelin, threshold_apop, threshold_stopmy)

reshaped_myelin = reshape(myelin.state, myelin.oligo_dim^2, myelin.oligo_counter);
indicate_destroyed_myelin = reshaped_myelin == 0;

total_damage_to_oligo = sum(indicate_destroyed_myelin); 

index1 = find(total_damage_to_oligo >= threshold_stopmy);

myelin.oligo_state(index1) = 2; % records oligo that can no longer myelinate 

index2 = find(total_damage_to_oligo >= threshold_apop); 

myelin.oligo_state(index2) = 0; % records oligo that is dead 

reshaped_myelin(:,index2) = 0;

myelin.state = reshape(reshaped_myelin,1,[]);

