function [myelin] = DEMYELINATION(myelin,C3)

% C3: create tall 2 column array where col 1 is x pos, col 2 is y pos 
C3_pos = [C3.x; C3.y]' ;

% Find where C3 cells occupy same place as myelin 
[~,index_my,~] = intersect(myelin.colarray,C3_pos,'rows'); 

% Degrade the health of occupied myelin by 1 grade 
myelin.state(index_my) = myelin.state(index_my) - 1; 
   
% Initialise the time since the myelin was degraded   
myelin.timer(index_my) = -1; 

% For efficiency, all occupied myelin is degraded, including state 0.  
% Set already dead myelin back to state '0' as that is the lowest state  
myelin.state(myelin.state < 0) = 0; 

end 