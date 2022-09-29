% this script should run on its own
close all; clear; clc

%--------------------------------------------------------------------------
% initial figure tasks
%--------------------------------------------------------------------------
commonFigureSetup;
commonFigureProperties;
hf = figure; % create a new figure and save handle
hf.Color = 'w'; % change the figure background color
hold on % do this once!

%--------------------------------------------------------------------------
% FIGURE-SPECIFIC CODE
%--------------------------------------------------------------------------
% data path
datapath = mfoldername(mfilename('fullpath'),'data');

% load data
data = load(fullfile(datapath,"figure2data"));

% plot the data
hp = surf(data.X,data.Y,data.Z,...
    "DisplayName","Peaks");

% change view
view(45,30)

% create axis labels
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')

% create legend (hl used in legend below)
hl = legend;
hl.Location = 'best';

% figure save name (used in export below)
savename = strcat("figure2");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'axes-z'; commonFigureTasks;
taskflag = 'legend'; commonFigureTasks;
taskflag = 'export'; commonFigureTasks;