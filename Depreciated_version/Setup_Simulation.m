 function  [funcs,PB,sizing,borders,myelin] = Setup_Simulation

%% Setup oligos
myelin.oligo_dim = 5; % dimension of oligo block 

if myelin.oligo_dim == 5 || myelin.oligo_dim == 10
     sizing.myelinwidth = 100; 
     sizing.domainheight = 300; 
elseif myelin.oligo_dim == 3
    sizing.myelinwidth = 99; 
    sizing.domainheight = 300; 
elseif myelin.oligo_dim == 7
        sizing.myelinwidth = 105; 
        sizing.domainheight = 301; 
end
%% Setup regions
% Define the size of each region 
sizing.cnswidth = 5; 
sizing.bloodwidth = 3; 
sizing.domainwidth = sizing.cnswidth + sizing.bloodwidth + sizing.myelinwidth; 

% Establish the borders of each region 
borders.leftblood = 0; 
borders.rightblood = borders.leftblood + sizing.bloodwidth; 
borders.leftcns = borders.rightblood; 
borders.rightcns = borders.leftcns + sizing.cnswidth; 
borders.leftmyelin = borders.rightcns; 
borders.rightmyelin = borders.leftmyelin + sizing.myelinwidth -1; 

%% Setup myelin
myelin.ybound = sizing.domainheight:-1:1; 
myelin.xbound =  borders.leftmyelin:1:borders.rightmyelin; 
[myelin.xmatrix,myelin.ymatrix] = meshgrid(myelin.xbound,myelin.ybound); 
myelin.xarray = reshape(myelin.xmatrix,1,[]);
myelin.yarray = reshape(myelin.ymatrix,1,[]);
myelin.allgone = 0;
myelin.depleted = 0; 

myelin.grades = 4; % this is n - 1 (code includes a zero state)
myelin.healtime = 25;
%% Setup PB
% PB properties
PB.leavingprob = 0.1; 
PB.enteringprob = 0; 
PB.rightside = borders.rightblood + 1;
PB.leftside = borders.rightblood;
   
%% Setup C3 biasing 
funcs.C3_aggression = 0.5 ; % C3 aggression + 0.25 = max bias towards myelin 
funcs.C3_bias = @(distance) funcs.C3_aggression*exp(-0.025*distance); 

