function [myelin, C1, C2, C3, tau] = Setup_Agents( ...
    myelin, borders, sizing)

% SETUP_AGENTS
%
% Initialises agent populations (C1, C2, C3) and constructs the myelin
% structure, including oligodendrocyte grouping and state variables.
%
% INPUTS:
%   myelin  - myelin structure from Setup_Simulation
%   borders - spatial region boundaries
%   sizing  - domain dimensions
%
% OUTPUTS:
%   myelin - updated myelin structure with oligodendrocyte mapping
%   C1     - Cell type 1 (peripheral blood agents)
%   C2     - Cell type 2 (CNS agents)
%   C3     - Cell type 3 (initially empty population)
%   tau    - timestep duration (minutes)

%% Time Parameters
tau = 20;                 % Timestep duration (minutes)
mins_per_day = 1440;

%% CELL TYPE 1 (C1)
C1.init      = 5; % initial number of cells 
C1.deathrate = 0.35 / (mins_per_day / tau);

C1.age = zeros(1, C1.init);
C1.x   = zeros(1, C1.init);
C1.y   = zeros(1, C1.init);

for ii = 1:C1.init % initialise spatial positions 
    C1.x(ii) = randi([borders.leftblood + 1, borders.rightblood - 1]);
    C1.y(ii) = randi([1, sizing.domainheight - 1]);
end

%% CELL TYPE 2 (C2)
C2.init = 45; % initial number of cells 

C2.x = zeros(1, C2.init);
C2.y = zeros(1, C2.init);

for ii = 1:C2.init % initialise spatial positions
    C2.x(ii) = randi([borders.leftcns + 1, borders.rightcns - 1]);
    C2.y(ii) = randi([1, sizing.domainheight - 1]);
end

%% CELL TYPE 3 (C3) 
C3.x = []; % set up empty arrays for when they are created during the simulation 
C3.y = [];
C3.age = [];

C3.deathrate = C1.deathrate;

%% MYELIN & OLIGODENDROCYTES
% Divide domain into oligodendrocyte-sized subregions

num_subheight = sizing.domainheight / myelin.oligo_dim;
num_subwidth  = sizing.myelinwidth  / myelin.oligo_dim;

C_x = mat2cell( ...
    myelin.xmatrix, ...
    repmat(myelin.oligo_dim, num_subheight, 1), ...
    repmat(myelin.oligo_dim, num_subwidth, 1));

C_y = mat2cell( ...
    myelin.ymatrix, ...
    repmat(myelin.oligo_dim, num_subheight, 1), ...
    repmat(myelin.oligo_dim, num_subwidth, 1));

% Initialise arrays
myelin.xarray        = [];
myelin.yarray        = [];
myelin.oligo_tracker = [];
myelin.oligo_counter = 0;

% Construct oligodendrocyte mapping
for ii = 1:num_subheight
    for jj = 1:num_subwidth
        
        myelin.oligo_counter = myelin.oligo_counter + 1;
        
        myelin.xarray = [myelin.xarray, ...
            reshape(C_x{ii, jj}, 1, [])];
        
        myelin.yarray = [myelin.yarray, ...
            reshape(C_y{ii, jj}, 1, [])];
        
        myelin.oligo_tracker = [ ...
            myelin.oligo_tracker, ...
            myelin.oligo_counter * ...
            ones(1, myelin.oligo_dim^2)];
    end
end

%% Myelin State Initialisation
myelin.numpieces        = length(myelin.xarray);

myelin.state            = myelin.grades * ...
                          ones(1, myelin.numpieces);

myelin.timer            = -1 * ones(1, myelin.numpieces);

myelin.fullyhealthyquant = ...
    myelin.grades * myelin.numpieces;

myelin.colarray = [myelin.xarray; myelin.yarray]';

% Oligodendrocyte states (1 = active)
myelin.oligo_state = ones(1, myelin.oligo_counter);

end