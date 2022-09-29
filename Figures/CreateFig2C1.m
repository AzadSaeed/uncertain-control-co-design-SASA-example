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
Pcase = 1;
foldername = strcat(solpath,'Case',string(Pcase));


%--------------------------------------------------------------------------
% FIGURE-SPECIFIC CODE
%--------------------------------------------------------------------------


Name = "sol_MC_N_Stc_MCS.mat";
load(strcat(foldername,filesep,Name));

% plot controls
hp1c = plot(sol.T,sol.U(:,:),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.green(3,:),'MarkerEdgeColor',C.green(3,:));


% Adjust the tick label
ax = gca;
ax.FontSize = fonttick;

% create axis labels
xlabel('$\textrm{Time}$','FontSize',fontlabel)
ylabel('$\textrm{States}$','FontSize',fontlabel)


% figure save name (used in export below)
savename = strcat(pwd,filesep,'Figures',filesep,"Fig2C1.pdf");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
exportgraphics(gca,savename);


