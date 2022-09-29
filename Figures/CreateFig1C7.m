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
GF = 0.1:0.1:5;

for i = 1:length(GF)

    filepath = erase(pwd,'Figures');
    solpath = strcat(filepath,'\Solutions\');
    Pcase = 7;
    foldername = strcat(solpath,'Case',string(Pcase));


    %--------------------------------------------------------------------------
    % FIGURE-SPECIFIC CODE
    %--------------------------------------------------------------------------


    Name = strcat("sol_MC_SH_Det_WCpoly",num2str(GF(i)),'.mat');
    load(strcat(foldername,filesep,Name));

    F(i) = sol.Y(end,1);
    k(i) = sol.P;
end


% plot state 1
fig = gca;
set(fig,'defaultAxesColorOrder',[C.orange(10,:); C.bluegrey(10,:)]);
yyaxis left
hp1c = plot(GF,F,'.',...
    'linewidth',linewidth,'markersize',16,...
    'Color',C.purple(10,:),'MarkerEdgeColor',C.purple(10,:));
ylabel('$\textrm{objective}$','FontSize',fontlabel)
ylim([-0.1, 0.35])
hold on

yyaxis right
hp1c = plot(GF,k,'.',...
    'linewidth',linewidth,'markersize',16,...
    'Color',C.orange(7,:),'MarkerEdgeColor',C.orange(7,:));
ylabel('$\mu_{k}$','FontSize',fontlabel)


% Adjust the tick label
ha = gca;
xlim([0.01 5])
ha.YAxis(1).Color = C.purple(10,:);
ha.YAxis(2).Color = C.orange(7,:);
ha.FontSize = fonttick;

% create axis labels
xlabel('$F$','FontSize',fontlabel)


% figure save name (used in export below)
savename = strcat(pwd,filesep,'Figures',filesep,"Fig1C7.pdf");

%--------------------------------------------------------------------------
% additional tasks
%--------------------------------------------------------------------------
ha.XAxis.Color = 'k'; % change the x axis color to black
ha.XAxis.FontSize = fonttick; % change x tick font size
ha.YAxis(1).FontSize = fonttick; % change y tick font size
ha.YAxis(2).FontSize = fonttick; % change y tick font size
ha.XAxis.Label.FontSize = fontlabel; % change x label font size
ha.YAxis(1).Label.FontSize = fontlabel; % change y label font size
ha.YAxis(2).Label.FontSize = fontlabel; % change y label font size
ha.Layer = 'top'; % place the axes on top of the data
ha.LineWidth = 1; % increase line width
taskflag = 'legend'; commonFigureTasks;
exportgraphics(gca,savename);


