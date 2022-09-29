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
% Directory settings
%--------------------------------------------------------------------------

filepath = erase(pwd,'Figures');
solpath = strcat(filepath,'\Solutions\');
Pcase = 3;
foldername = strcat(solpath,'Case',string(Pcase));


%--------------------------------------------------------------------------
% FIGURE-SPECIFIC CODE
%--------------------------------------------------------------------------


Name = "sol_SC_SH_Det_WC.mat";
load(strcat(foldername,filesep,Name));

% plot state 1
hp1c = plot(sol.T,sol.Y(:,1),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.red(3,:),'MarkerEdgeColor',C.red(3,:));

% plot state 2
hp1c = plot(sol.T,sol.Y(:,2),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.blue(3,:),'MarkerEdgeColor',C.blue(3,:));


% Adjust the tick label
ax = gca;
ax.FontSize = fonttick;

% create axis labels
xlabel('$\textrm{Time}$','FontSize',fontlabel)
ylabel('$\textrm{States}$','FontSize',fontlabel)


% figure save name (used in export below)
savename = strcat(pwd,filesep,'Figures',filesep,"Fig1C3.pdf");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
exportgraphics(gca,savename);


