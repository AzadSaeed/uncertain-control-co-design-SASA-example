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
% create data
[X,Y,Z] = peaks(250);

% plot the data
hp = contourf(X,Y,Z,50,...
    'DisplayName','Peaks');

% create axis labels
xlabel('$x$')
ylabel('$y$')

% create legend (hl used in legend below)
hl = legend;

% figure save name (used in export below)
savename = "figure3";

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
taskflag = 'export'; commonFigureTasks;