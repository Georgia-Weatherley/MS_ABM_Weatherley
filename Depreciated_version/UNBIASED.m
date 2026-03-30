function [xpos,ypos] = UNBIASED(xpos,ypos,PB)

delta = 1; % step size 
move_probs = [0.25 0.25 0.25 0.25]; % movement probabilities
cumulative_probs = cumsum(move_probs);
num_cells = length(xpos);

r_C = rand(1,num_cells); % generate random number for each cell

west_mask = r_C < cumulative_probs(1); % mask of cells to move west
west_mask(xpos == 1) = 0; 
if sum(west_mask) > 0 
    r_BBB_west = rand(1,num_cells);
    accept_mask = r_BBB_west <= PB.enteringprob; % see if cell is accepted to move across BBB
    west_mask(west_mask == 1 & xpos == PB.rightside & accept_mask ~= 1) = 0;
    xpos(west_mask)  = xpos(west_mask)  - delta; % move those cells west by delta 
end
 
east_mask = r_C < cumulative_probs(2) & r_C >= cumulative_probs(1) ; % mask of cells to move east 
if sum(east_mask) > 0 
    r_BBB_east = rand(1,num_cells);
    accept_mask = r_BBB_east <= PB.leavingprob; % see if cell is accepted to move across BBB
    east_mask(east_mask == 1 & xpos == PB.leftside & accept_mask ~= 1) = 0;
    xpos(east_mask)  = xpos(east_mask)  + delta; % move those cells east by delta 
end

south_mask = r_C < cumulative_probs(3) & r_C >= cumulative_probs(2); % mask of cells to move south 
ypos(south_mask) = ypos(south_mask) - delta; % move those cells south by delta 

north_mask = r_C >= cumulative_probs(3) ; % mask of cells to move north
ypos(north_mask) = ypos(north_mask) + delta; % move those cells north by delta 

end 
     