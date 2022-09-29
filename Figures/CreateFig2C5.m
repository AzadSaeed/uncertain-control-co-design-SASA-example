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

for i = 1:8

    filepath = erase(pwd,'Figures');
    solpath = strcat(filepath,'\Solutions\');
    Pcase = 5;
    foldername = strcat(solpath,'Case',string(Pcase));


    %--------------------------------------------------------------------------
    % FIGURE-SPECIFIC CODE
    %--------------------------------------------------------------------------


    Name = strcat("sol_MC_SH_Det_WCpoly",num2str(i),'.mat');
    load(strcat(foldername,filesep,Name));

    % plot state 1
    hp1c = plot(sol.T,sol.U(:,1),'.',...
        'linewidth',linewidth,'markersize',12,...
        'Color',C.green(3,:),'MarkerEdgeColor',C.green(3,:));

end

% Adjust the tick label
ax = gca;
ax.FontSize = fonttick;

% create axis labels
xlabel('$\textrm{Time}$','FontSize',fontlabel)
ylabel('$\textrm{States}$','FontSize',fontlabel)


% figure save name (used in export below)
savename = strcat(pwd,filesep,'Figures',filesep,"Fig2C5.pdf");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
taskflag = 'axes'; commonFigureTasks; %#ok<*NASGU>
taskflag = 'legend'; commonFigureTasks;
exportgraphics(gca,savename);


