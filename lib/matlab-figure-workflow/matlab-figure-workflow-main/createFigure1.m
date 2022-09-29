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
x = linspace(-2*pi,2*pi,50);
y = exp(-0.3*x).*sin(x);

% plot the data
hp = plot(x,y,'.-',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.red(9,:),'MarkerEdgeColor',C.blue(4,:),...
    'DisplayName','$y = e^{-3x/10} \sin(x)$');

% create axis labels
xlabel('$x$')
ylabel('$y$')

% create legend (hl used in legend below)
hl = legend;

% figure save name (used in export below)
savename = strcat("figure1","-n",string(length(x)));

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
taskflag = 'export'; commonFigureTasks;