function [funcs, PB, sizing, borders, myelin] = Setup_Simulation

% SETUP_SIMULATION
%
% Initialises the spatial domain, region boundaries, myelin grid,
% and model parameters for the simulation.
%
% OUTPUTS:
%   funcs   - function handles and behavioural parameters
%   PB      - probabilistic boundary parameters
%   sizing  - domain and region dimensions
%   borders - spatial region boundaries
%   myelin  - myelin grid and associated properties


%% Oligodendrocyte Setup

myelin.oligo_dim = 5;   % Dimension of oligo block

% Domain sizing based on oligo dimension
if myelin.oligo_dim == 5 || myelin.oligo_dim == 10
    sizing.myelinwidth  = 100;
    sizing.domainheight = 300;

elseif myelin.oligo_dim == 3
    sizing.myelinwidth  = 99;
    sizing.domainheight = 300;

elseif myelin.oligo_dim == 7
    sizing.myelinwidth  = 105;
    sizing.domainheight = 301;
end

%% Region Setup

% Define widths of each region
sizing.cnswidth   = 5;
sizing.bloodwidth = 3;

% Total domain width
sizing.domainwidth = sizing.cnswidth + ...
                     sizing.bloodwidth + ...
                     sizing.myelinwidth;

% Establish region boundaries
borders.leftblood   = 0;
borders.rightblood  = borders.leftblood + sizing.bloodwidth;

borders.leftcns     = borders.rightblood;
borders.rightcns    = borders.leftcns + sizing.cnswidth;

borders.leftmyelin  = borders.rightcns;
borders.rightmyelin = borders.leftmyelin + sizing.myelinwidth - 1;

%% Myelin Grid Setup
myelin.ybound = sizing.domainheight:-1:1;
myelin.xbound = borders.leftmyelin:1:borders.rightmyelin;

[myelin.xmatrix, myelin.ymatrix] = meshgrid( ...
    myelin.xbound, myelin.ybound);

myelin.xarray = reshape(myelin.xmatrix, 1, []);
myelin.yarray = reshape(myelin.ymatrix, 1, []);

% Myelin state indicators
myelin.allgone   = 0;
myelin.depleted  = 0;

% Myelin properties
myelin.grades   = 4;   % n - 1 (includes zero state)
myelin.healtime = 25;

%% Probabilistic Boundary (PB)
PB.leavingprob  = 0.1;
PB.enteringprob = 0;

PB.rightside = borders.rightblood + 1;
PB.leftside  = borders.rightblood;

%% C3 Biasing Functions
funcs.C3_aggression = 0.5;   % Bias parameter

% Bias function: increases attraction toward myelin based on distance
funcs.C3_bias = @(distance) ...
    funcs.C3_aggression * exp(-0.025 * distance);

end