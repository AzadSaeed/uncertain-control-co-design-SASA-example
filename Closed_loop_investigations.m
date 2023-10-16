clear; clc
close all


% Plot histogram for final state for all cases
for i = 1:1
    [Mu(i), Sigma(i)] = CreateHis(i)
end


% Report Settling time
c = 0.02; % criterial for settling time (2%)
for i = 1:7
    T(i) = FindSettlingTime(i,c);
end


Dir = strcat(pwd,filesep,'Solutions',filesep, 'Case8');
if ~isfolder(Dir)
    mkdir(Dir);
end



% Define the cases to study, with some options
TestCases  = {'Det', 'MCS', 'WCR'};
opts4      = pidtuneOptions('DesignFocus','referencetracking','PhaseMargin', 90);
opts       = odeset('RelTol',1e-11,'AbsTol',1e-12);
paropts    = parforOptions(gcp, 'MaxNumWorkers', Inf);



% Deterministic Problem Data
J          = 1;                       % Inertial ratioe
x0_nom     = [0; 0];                  % Initial state in nominal case
t_end      = 5;                       % final time
t_span     = [0 t_end];               
tt         = linspace(0,t_end,1000)';  % time 
N          = 5;                       % Number of samples in MCS
P          = [10, 50, 90, 100];       % Percentiles for MCS






for idx = 3:length(TestCases)


    % Load data associated with each case in TestCases
    [t,u,y,K] = LoadData(TestCases{idx});



    % Create the reference
    % Ref{1,1} ..... relative displacement
    % Ref{1,2} ..... relative velocity
    % Ref{1,3} ..... nominal control
    Ref = CreateRef(t,u,y,TestCases{idx},P);
  
    C_nom.Kp = 0;              % Gain for tracking displcament trajectory
    C_nom.Kd = 10^6;           % Gain for tracking velocity trajectory

    % Solve with uncertainties
    if idx ==2

        for i = 1:size(Ref,1)
            Ref_            = {Ref{i,1}, Ref{i,2}, Ref{i,3}};
            [Y1_u, Y2_u, U] = solve_U_odes(K,J,x0_nom,C_nom, Ref_, N, paropts,opts, t_span,idx,tt);
            Y_u             =[Y1_u,Y2_u,U];
            [R1,R2]         = ExtendRef(tt,Ref);
            RefVal          = {R1,R2};

            % Plots
            P = [10, 50, 90, 100];
            Name = strcat('MCS_',num2str(P(i)));
            varclass = {'State', 'Control'};
            createplots(tt,[],Y_u,RefVal,strcat(Dir,filesep,Name),varclass{1})
            createplots(tt,[],U,RefVal,strcat(Dir,filesep,Name),varclass{2})

            save(strcat(Dir,filesep,'Y1_u_',Name), 'Y1_u')
            save(strcat(Dir,filesep,'Y2_u_',Name), 'Y2_u')
            save(strcat(Dir,filesep,'U_', Name), 'U')

        end




    else

        [Y1_u, Y2_u, U] = solve_U_odes(K,J,x0_nom,C_nom, Ref, N, paropts,opts, t_span,idx,tt);
        Y_u             = [Y1_u,Y2_u,U];
        [R1,R2]         = ExtendRef(tt,Ref);
        RefVal          = {R1,R2};

        % Plots
        Name = {'Det','MCS', 'WCR'};
        varclass = {'State', 'Control'};
        createplots(tt,[],Y_u,RefVal,strcat(Dir,filesep,Name{idx}),varclass{1})
        createplots(tt,[],U,RefVal,strcat(Dir,filesep,Name{idx}),varclass{2})

        save(strcat(Dir,filesep,'Y1_u_',Name{idx}), 'Y1_u')
        save(strcat(Dir,filesep,'Y2_u_',Name{idx}), 'Y2_u')
        save(strcat(Dir,filesep,'U_', Name{idx}), 'U')

    end


close all;
end







%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function varargout = LoadData(TestCase)

switch TestCase

    case 'Det'

        Det = load(strcat(pwd,filesep,'Solutions',filesep,'Case0',filesep,'sol_SC_N_Det_Det'));
        t   = Det.sol.T;
        u   = Det.sol.U;
        y   = Det.sol.Y;
        k   = Det.sol.xopt;
        varargout = {t,u,y,k};

    case 'MCS'

        MCS = load(strcat(pwd,filesep,'Solutions',filesep,'Case1',filesep,'sol_MC_N_Stc_MCS'));
        t   = MCS.sol.T;
        u   = MCS.sol.U;
        y1  = MCS.sol.Y1;
        y2  = MCS.sol.Y2;
        y   = [y1,y2];
        k   = MCS.sol.xopt;
        varargout = {t,u,y,k};

    case 'WCR'

        WCR = load(strcat(pwd,filesep,'Solutions',filesep,'Case5',filesep,'sol_MC_N_Det_WCpoly7'));       
        WCR.sol.k = WCR.sol.xopt + 0.6;

        t = WCR.sol.T;
        u = WCR.sol.U;
        y = WCR.sol.Y;
        k = WCR.sol.k;
        varargout = {t,u,y,k};

end

end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function Ref =  CreateRef(t,u,y,test_case,P)


switch test_case

    case 'Det' 

        Ref_y1 = griddedInterpolant(t,y(:,1),'spline');
        Ref_y2 = griddedInterpolant(t,y(:,2),'spline');
        Ref_u  = griddedInterpolant(t,u,'linear');
        Ref = {Ref_y1,Ref_y2,Ref_u};

    case 'MCS'

        n   = size(u,2);
        y1  = y(:,1:n);
        y2  = y(:,n+1:2*n);

        % Create multiple references based on various percentiles
        Ref = cell(length(P),3);

        for idx = 1:length(P)

            y1_    = prctile(y1',P(idx));
            y2_    = prctile(y2',P(idx));
            u_     = prctile(u',P(idx));
            Ref_y1 = griddedInterpolant(t,y1_,'spline');
            Ref_y2 = griddedInterpolant(t,y2_,'spline');
            Ref_u  = griddedInterpolant(t,u_,'linear');
            Ref(idx,:)    = {Ref_y1, Ref_y2, Ref_u};

        end


    case 'WCR'

        Ref_y1 = griddedInterpolant(t,y(:,1),'spline');
        Ref_y2 = griddedInterpolant(t,y(:,2),'spline');
        Ref_u  = griddedInterpolant(t,u,'linear');
        Ref    = {Ref_y1,Ref_y2,Ref_u};
end

end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function [A,B,C,D] = DynModel(k,J)

A = [0, 1; -k/J, 0];
B = [0; 1/J];
C = [1,1];
D = 0;


end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function dx  = sys(t,x, K, J, Gains,R)

Kp = Gains.Kp;
Kd = Gains.Kd;

[A,B,~,~] = DynModel(K,J);



%% Using relative displacement and velocity as reference;
if t <=1
    rx   = R{1,1}(t);           % reference trajectory - relative displacement
    rv   = R{1,2}(t);           % reference trajectory - relative velocity
    U    = R{1,3}(t);           % nominal control trajectory
elseif t>1
    rx   = R{1,1}(1);           % For simulation more than 1 s - keep velocity at 0
    rv   = R{1,2}(1);           % keep position at its final value
    U    = K*x(1);              % Calculate control based on current position
end

ex   = -(x(1) - rx);            % error in relative displacement
ev   = -(x(2) - rv);            % error inn relative velocity


u = Kp*ex + Kd*ev + U;


% Check for control saturation
umax = 1;
if u>umax
    u = umax;
elseif u<-umax
    u = -umax;
end

dx = A*x + B*u;


end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function U = CalcU(t,x, Gains, k, R)

Kp   = Gains.Kp;
Kd   = Gains.Kd;


if t(end)>1

    idx_          = find(t>=1);
    Idx           = idx_(1);
    rx            = R{1,1}(t);
    rv            = R{1,2}(t);
    U             = R{1,3}(t);
    rx(Idx:end)   = R{1,1}(1);
    rv(Idx:end)   = R{1,2}(1);
    U(Idx:end)    = k*x(Idx:end,1);

else

    rx            = R{1,1}(t);
    rv            = R{1,2}(t);
    U             = R{1,3}(t);

end



ex   = -(x(:,1) - rx);            % error in relative displacement
ev   = -(x(:,2) - rv);            % error inn relative velocity


u    = Kp*ex + Kd*ev + U;


umax = 1;
ro_u    = u>umax;
u(ro_u) = umax;
ro_l    = find(u<-umax);
u(ro_l) = -umax;

U = u;


end




%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function detplots(tt,Y,Ref,filename)


ha = figure; hold on; rgb_color = materialColors;

% Potentially re-create references for t-end > 1
[R1,R2] = ExtendRef(tt,Ref);

% PLot reference
plot(tt, R1,'LineWidth',3, 'Color',rgb_color.grey(7,:))
plot(tt, R2,'LineWidth',3, 'Color',rgb_color.grey(7,:))

% PLot states
plot(tt, Y{1,1}(tt),'LineWidth',3, 'Color',rgb_color.red(7,:))
plot(tt, Y{1,2}(tt),'LineWidth',3, 'Color',rgb_color.blue(7,:))



% Plot Control
DATA_Color  = [rgb_color.blue(6,:);rgb_color.orange(6,:);rgb_color.purple(6,:)...
                ;rgb_color.teal(6,:);rgb_color.lime(6,:);rgb_color.green(6,:)];

plot(tt, Y{1,3}(tt),'LineWidth',3, 'Color',DATA_Color(2,:))

xlabel('time')
ylabel('states and controls')

exportgraphics(ha,strcat(filename,'.pdf'),'Resolution',400)


end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



function [Y1_u, Y2_u, U] = solve_U_odes(K,J,x0_nom,C_nom, Ref, N, paropts,opts, t_span,idx, tt)


% Uncertain problem data
SK         = 0.2;    % Standard Deviation for Plant Variable (k)
SJ         = 0.15;   % Standard Deviation for Fixed Problem Parameter (J)
SY20       = 0.03;   % Standard Deviation for Initial Boundary Condition (Y0)
k_idx      = 3;


% Uncertain samples

Mu_K   = K-0.6;
Mu_J   = J;
Mu_Y02 = x0_nom(2);

if idx ==1 || idx ==2

    rng(0.854)
    K_smpl = randn(N,1)*SK + Mu_K;
    rng(0.654)
    J_smpl = randn(N,1)*SJ + Mu_J;
    rng(0.128)
    Y_smpl = randn(N,1)*SY20 + Mu_Y02;

elseif idx ==3

    % samples should be created such that Mu_k obtained from UCCD solution
    % is the mean value. Accordingly, since we are using the solution at
    % vertex 7, the interval for stiffness is [5.8-0.6 5.8+0.6] = [5.2 6.4]
    k_l = Mu_K - k_idx*SK;
    k_h = Mu_K + k_idx*SK;
    rng(0.854)
    K_smpl = k_l + (k_h -k_l).*rand(N,1);

    J_l = Mu_J - k_idx*SJ;
    J_h = Mu_J + k_idx*SJ;
    rng(0.654)
    J_smpl = J_l + (J_h -J_l).*rand(N,1);

    Y_l = -k_idx*SY20;
    Y_h = -k_idx*SY20;
    rng(0.128)
    Y_smpl = Y_l + (Y_h -Y_l).*rand(N,1);


end


Y1_u    = zeros(length(tt),N);
Y2_u    = zeros(length(tt),N);
x0_nom1 = x0_nom(1);
U       = zeros(length(tt),N);


% Run simulations
parfor (i = 1:N,paropts)
%for i = 1:N

    [t_u, y_u]  = ode15s(@(t,x)sys(t,x,K_smpl(i),J_smpl(i),C_nom,Ref),...
        t_span, [x0_nom1;Y_smpl(i)], opts);

    F   = griddedInterpolant(t_u,y_u(:,1),'spline');
    Y1_u(:,i) = F(tt);
    F   = griddedInterpolant(t_u,y_u(:,2),'spline');
    Y2_u(:,i) = F(tt);


end


% Calculate U 
parfor (i = 1:N,paropts)
    U(:,i)           = CalcU(tt,[Y1_u(:,i),Y2_u(:,i)], C_nom, K_smpl(i), Ref);
end


%%%%%%   Additional investigaitons %%%%%%
% kset = [-k_idx*SK, k_idx*SK];
% Jset = [-k_idx*SJ, k_idx*SJ];
% SY20set = [-k_idx*SY20, k_idx*SY20];
% ll= 1;
% for ii=1:2
%     for jj=1:2
%         for kk=1:2
%             vert(ll,:) =[kset(ii),Jset(jj),SY20set(kk)]; 
%             ll = ll + 1;
%         end
%     end
% end

% for hh = 1:8
%     [t_u, y_u]  = ode15s(@(t,x)sys(t,x,Mu_K+vert(hh,1),Mu_J+vert(hh,2),C_nom,Ref),...
%         t_span, [x0_nom1;vert(hh,3)], opts);
%     F   = griddedInterpolant(t_u,y_u(:,1),'spline');
%     Y1_u_v(:,hh) = F(tt);
%     F   = griddedInterpolant(t_u,y_u(:,2),'spline');
%     Y2_u_v(:,hh) = F(tt);
%     clear t_u y_u
% end

% for hh = 1:8
%     U_v(:,hh)           = CalcU(tt,[Y1_u_v(:,hh),Y2_u_v(:,hh)], C_nom, vert(hh,1), Ref);
% end


% [R1,R2] = ExtendRef(tt,Ref);
% rgb_color = materialColors;

% for i=1:N
%     figure(1) 
%     hold on
%     plot(tt,U(:,i),'LineWidth',3,'Color',rgb_color.green(5,:))
% end


% for i=1:N
%     figure(2)
%     hold on
%     plot(tt,Y1_u(:,i),'LineWidth',3,'Color',rgb_color.red(3,:))
% end
% for i=1:8
%     plot(tt,Y1_u_v(:,i),'LineWidth',3,'Color',rgb_color.grey(10,:))
% end
% 
% 
% for i=1:N
%     figure(3)
%     plot(tt,Y2_u(:,i),'LineWidth',3,'Color',rgb_color.blue(3,:))
%     hold on
% end
% for i=1:8
%     plot(tt,Y1_u_v(:,i),'LineWidth',3,'Color',rgb_color.grey(10,:))
% end
% plot(tt,R2,'LineWidth',3,'Color',rgb_color.grey(5,:))


end




%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [R1, R2] = ExtendRef(t,Ref)

if t(end)>1

    idx_          = find(t>=1);
    Idx           = idx_(1);
    R1            = Ref{1,1}(t);
    R2            = Ref{1,2}(t);
    R1(Idx:end)   = Ref{1,1}(1);
    R2(Idx:end)   = Ref{1,2}(1);

else

    R1  = Ref{1,1}(tt);
    R2  = Ref{1,2}(tt);

end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function createplots(t,Y_nom,Y_u,Refval,Name, varclass)

n     = size(Y_u,2)/3;
Y1    = Y_u(:,1:n);
Y2    = Y_u(:,n+1:2*n);
U     = Y_u(:,2*n+1:3*n);



% settings for plots;
set(0,'DefaultTextInterpreter','latex'); % change text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change tick interpreter
rgb_color = materialColors;
hf  = figure; hf.Color = 'W'; hold on


switch varclass

    case 'State'

        % Plot uncertain responce of Y1
        Y1_min = min(Y1,[],2);
        Y1_max = max(Y1,[],2);
        f1 = patch([t;flipud(t)],[Y1_min; flipud(Y1_max)], rgb_color.orange(7,:),...
            'EdgeColor', rgb_color.orange(7,:),'FaceAlpha',0.4,'LineStyle', 'none');

        p1 = plot(t,prctile(Y1',10),'LineWidth',3, 'Color',rgb_color.grey(5,:));

        p2 = plot(t,prctile(Y1',50),'LineWidth',3, 'Color',rgb_color.grey(6,:));

        p3 = plot(t,prctile(Y1',90),'LineWidth',3, 'Color',rgb_color.grey(7,:));


        % Plot uncertain responce of Y2
        Y2_min = min(Y2,[],2);
        Y2_max = max(Y2,[],2);
        f2 = patch([t;flipud(t)],[Y2_min; flipud(Y2_max)], rgb_color.orange(7,:),...
            'EdgeColor', rgb_color.orange(7,:),'FaceAlpha',0.4,'LineStyle', 'none');

        p1 = plot(t,prctile(Y2',10),'LineWidth',3, 'Color',rgb_color.grey(5,:));

        p2 = plot(t,prctile(Y2',50),'LineWidth',3, 'Color',rgb_color.grey(6,:));

        p3 = plot(t,prctile(Y2',90),'LineWidth',3, 'Color',rgb_color.grey(7,:));


        % Plot reference trajectories
 
        f3 = plot(t,Refval{1,1},'LineWidth',3, 'Color',rgb_color.red(7,:));
        f4 = plot(t,Refval{1,2},'LineWidth',3, 'Color',rgb_color.blue(7,:));


        ylabel('Relative displacement (m)')
        %         legend([f1,f3],sprintf('Uncertain response to reference tracking based on %s CCD solution',Name),...
        %             'Reference trajectory', 'Location', 'northoutside' )
        % xlim([-0.12 0.7])

    case 'Control'

        for i = 1:size(U,2)
            % f1 = plot(t,Y_u(:,i),'.','MarkerSize',7, 'MarkerEdgeColor',...
            %      rgb_color.orange(7,:),'MarkerFaceColor', rgb_color.orange(7,:));
            f1 = scatter(t,Y_u(:,i),'o','filled','MarkerFaceColor', rgb_color.orange(7,:),...
                'MarkerEdgeColor', rgb_color.orange(7,:),'MarkerFaceAlpha',...
                0.2, 'MarkerEdgeAlpha',0.2);
        end


        ylabel('Control Effort ')
        %legend(f1,sprintf('Uncertain control effort based on %s CCD solution',Name), 'Location', 'northoutside' )
        ylim([-1.1, 1.1])

end

ha           = gca;
ha.YColor    = 'k';
ha.Layer     = 'top'; % place the axes on top of the data
ha.LineWidth = 2; % increase axis line width
xlabel('Time (sec)')




filename = strcat(Name,'_',varclass);
exportgraphics(ha,strcat(filename,'.pdf'),'Resolution',400)


end


%------------------------------------$-------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



function [Mu, Sigma] = CreateHis(i)

% settings for plots;
set(0,'DefaultTextInterpreter','latex'); % change text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change tick interpreter

rgb_color   = materialColors;
hf          = figure; hf.Color = 'W'; hold on

% Load solutions
name        = {'Y1_u_Det_P','Y1_u_Det_B','Y1_u_MCS_10','Y1_u_MCS_50','Y1_u_MCS_90','Y1_u_MCS_100','Y1_u_WCR'};
datapath    = strcat(pwd,filesep,'Solutions', filesep,'Case8',filesep, name{1,i});
Data         = load(datapath);


Edges       = linspace(0,0.5,50);
DATA_Color  = [rgb_color.blue(3,:);rgb_color.lightblue(5,:); rgb_color.orange(4,:);rgb_color.purple(4,:)...
                ;rgb_color.teal(4,:);rgb_color.lime(6,:);rgb_color.green(6,:)];


[filepath,name,ext] = fileparts(datapath) 
filename    = strcat(filepath,filesep,'obj',string(i));
alpha       = 0.7;


if i ==1 || i == 2 || i==7
    histogram(Data.Y1_u(end,:), Edges, 'FaceColor', DATA_Color(i,:),...
        'EdgeColor',DATA_Color(i,:), 'FaceAlpha', alpha,'EdgeAlpha',alpha,...
        'Normalization','probability');
    plot(mean(Data.Y1_u(end,:))*ones(1,10),linspace(0,1,10),'LineWidth',3, 'Color',rgb_color.grey(5,:))
    Mu_ = mean(Data.Y1_u(end,:));
    Sigma_ = std(Data.Y1_u(end,:));
else
    histogram(Data.Y1_u.Y1_u(end,:), Edges, 'FaceColor', DATA_Color(i,:),...
        'EdgeColor',DATA_Color(i,:), 'FaceAlpha', alpha,'EdgeAlpha',alpha,...
        'Normalization','probability');
    plot(mean(Data.Y1_u.Y1_u(end,:))*ones(1,10),linspace(0,1,10),'LineWidth',3, 'Color',rgb_color.grey(5,:))
    Mu_ = mean(Data.Y1_u.Y1_u(end,:));
    Sigma_ = std(Data.Y1_u.Y1_u(end,:));
end

xlabel('objective function value')
ylabel('probability')
xlim([0, 0.5])
ylim([0 1])
ha               = gca;
ha.YColor        = 'k';
ha.Layer         = 'top'; % place the axes on top of the data
ha.LineWidth     = 2; % increase axis line width
figHandle        = gca;
exportgraphics(figHandle,strcat(filename,'.pdf'),'Resolution',400)

Mu = Mu_;
Sigma = Sigma_;

end



function T = FindSettlingTime(i,c)

name        = {'Y2_u_Det_P','Y2_u_Det_B','Y2_u_MCS_10','Y2_u_MCS_50','Y2_u_MCS_90','Y2_u_MCS_100','Y2_u_WCR'};
datapath    = strcat(pwd,filesep,'Solutions', filesep,'Case8',filesep, name{1,i});
Data        = load(datapath);
tt          = linspace(0,5,1000)';  % time 

if i ==1 || i ==2 || i == 7
    Y2 = Data.Y2_u;
else
    Y2 = Data.Y2_u;
end


idx_ = 1000:-1:201; % time indices greater than 1
v_mean = zeros(length(idx_),1);
v_min = zeros(length(idx_),1);
v_max = zeros(length(idx_),1);

for kk = 1:length(idx_) 
    v_mean(kk,1) = mean(Y2(idx_(kk),:));
    v_min(kk,1) = min(Y2(idx_(kk),:));
    v_max(kk,1) = max(Y2(idx_(kk),:));
end

% If greater than 0.02%, then the system has not settled yet
% choose the first instance where the mean velocity is less than 0.02
L_min = find(v_min <= -c);
L_max = find(v_max >= c);
L_mean = find(v_mean>=c | v_mean<-c);

if isempty(L_min)
    T_(1) = 1;
else 
    IDX = L_min(1);
    t_idx = idx_(IDX);
    T_(1) = tt(t_idx);
end

if isempty(L_max)
    T_(2) = 1;
else
    IDX = L_max(1);
    t_idx = idx_(IDX);
    T_(2) = tt(t_idx);
end

if isempty(L_mean)
    T_(3) = 1;
else
    IDX = L_mean(1);
    t_idx = idx_(IDX);
    T_(3) = tt(t_idx);
end

T = max(T_);

plot(tt(201:end),flip(v_min),'.r', tt(201:end), flip(v_max),'.b', tt(201:end), flip(v_mean))
hold on
plot(tt(201:end), c*ones(800,1),'.k',tt(201:end), -c*ones(800,1),'.k')
close all
end

