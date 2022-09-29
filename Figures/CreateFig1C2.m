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
Pcase = 2;
foldername = strcat(solpath,'Case',string(Pcase));


%--------------------------------------------------------------------------
% FIGURE-SPECIFIC CODE
%--------------------------------------------------------------------------


Name = "sol_MC_N_Stc_gpc.mat";
load(strcat(foldername,filesep,Name));

% plot state 1
hp1c = plot(sol.T,sol.EY1(:,1),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.red(3,:),'MarkerEdgeColor',C.red(3,:));

% plot state 2
hp1c = plot(sol.T,sol.EY2(:,1),'.',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.blue(3,:),'MarkerEdgeColor',C.blue(3,:));

% Plot EY1 + Sigma Y1
Y1_l = sol.EY1 - sqrt(sol.VarY1);
Y1_u = sol.EY1 + sqrt(sol.VarY1);
p1 = patch([sol.T;flipud(sol.T)],[Y1_l;flipud(Y1_u)],...
    C.orange(7,:),'EdgeColor', C.orange(7,:),'FaceAlpha', 0.3,...
    'LineStyle','none');


% Plot EY2 + Sigma Y2
Y2_l = sol.EY2 - sqrt(sol.VarY2);
Y2_u = sol.EY2 + sqrt(sol.VarY2);
p2 = patch([sol.T;flipud(sol.T)],[Y2_l;flipud(Y2_u)],...
    C.orange(7,:),'EdgeColor', C.orange(7,:),'FaceAlpha', 0.3,...
    'LineStyle','none');


% Adjust the tick label
ax = gca;
ax.FontSize = fonttick;

% create axis labels
xlabel('$\textrm{Time}$','FontSize',fontlabel)
ylabel('$\textrm{States}$','FontSize',fontlabel)


% figure save name (used in export below)
savename = strcat(pwd,filesep,'Figures',filesep,"Fig1C2.pdf");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
exportgraphics(gca,savename);


