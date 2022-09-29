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

% plot state 1
hp1c = plot(sol.T,mean(sol.Y1(:,:),2),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.red(3,:),'MarkerEdgeColor',C.red(3,:));

% plot state 2
hp1c = plot(sol.T,mean(sol.Y2(:,:),2),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.blue(3,:),'MarkerEdgeColor',C.blue(3,:));


% plot pacth for y1
y1min = min(sol.Y1,[],2); y1max = max(sol.Y1,[],2);
p1 = patch([sol.T;flipud(sol.T)], [y1min; flipud(y1max)],...
    C.orange(7,:),'EdgeColor', C.orange(7,:), 'FaceAlpha', 0.3,'LineStyle','none');

% plot patch for y2
y2min = min(sol.Y2,[],2); y2max = max(sol.Y2,[],2);
p2 = patch([sol.T;flipud(sol.T)], [y2min; flipud(y2max)],...
    C.orange(7,:),'EdgeColor', C.orange(7,:), 'FaceAlpha', 0.3,'LineStyle','none');

% Adjust the tick label
ax = gca;
ax.FontSize = fonttick;

% create axis labels
xlabel('$\textrm{Time}$','FontSize',fontlabel)
ylabel('$\textrm{States}$','FontSize',fontlabel)


% figure save name (used in export below)
savename = strcat(pwd,filesep,'Figures',filesep,"Fig1C1.pdf");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
exportgraphics(gca,savename);


