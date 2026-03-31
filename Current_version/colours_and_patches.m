function [col,patching,plotstyle] = colours_and_patches(borders,sizing)

% COLOURS_AND_PATCHES 
%
% This function is used to define custom colours and consistent plot
% styles. It takes borders and sizing to define the regions of the domain
% as patches for spatial plotting 

%% Custom colours 
col.BBB = (1/255)*[244 205 204]; % Peripheral Blood 
col.PB = (1/255)*[244 132 132]; % BBB
col.CNS =  (1/255)*[208,209,230]; % colour of CNS

col.G4 = (1/255)*[4,90,141]; % Myelin grade
col.G3 = (1/255)*[5,112,176]; % Myelin grade
col.G2 = (1/255)*[54,144,192]; % Myelin grade
col.G1 = (1/255)*[116,169,207]; % Myelin grade
col.G0 = (1/255)*[166,189,219];%[212,185,218]; % Myelin grade

col.C1 = 1/255*[254,196,79] ; % C1 cells 
col.C3 = 1/255*[85,153,255]; %C2 cells 
col.C3edge = (1/255)*[233, 116, 81] ;
col.C2 = 1/255*[127,42,255]; %C3 cells 

%% Plot styles 
plotstyle.width = 4.5; 
plotstyle.headingsize = 22;
plotstyle.labelsize = 20;
plotstyle.font = 'Times';
plotstyle.cell_markersize = 4;
plotstyle.myelin_markersize = 3; 

%% Patches for regions of domain
patching.y = [0 0 sizing.domainheight sizing.domainheight]; 
patching.CNS_x = [borders.leftcns+0.5 borders.rightcns borders.rightcns borders.leftcns+0.5];
patching.CNS_y = [0 0 sizing.domainheight sizing.domainheight]; 
patching.PB_x = [borders.leftblood borders.rightblood+0.5 borders.rightblood+0.5 borders.leftblood]; 
patching.PB_y = [0 0 sizing.domainheight sizing.domainheight]; 
patching.myelin_x = [borders.rightcns sizing.domainwidth sizing.domainwidth borders.rightcns];
patching.myelin_y = [0 0 sizing.domainheight sizing.domainheight]; 