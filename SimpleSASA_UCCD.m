function varargout = SimpleSASA_UCCD(varargin)

% Set general options
ex_opts =@SimpleSASA_UCCD_Opts;
[auxdata, opts] = DTQP_standardizedinputs(ex_opts, []);

% Assign inputs
if nargin == 1
    pcase = varargin{1};           % case-study number
elseif nargin == 2
    pcase = varargin{1};           % case-study number
    auxdata.Pen = varargin{2};     % penalty vlaue 
elseif nargin == 3
    pcase = varargin{1};           % case-study number
    auxdata.Pen = varargin{2};     % penalty value
    auxdata.GF = varargin{3};      % geometric (size) factor for std
elseif nargin == 4
    pcase = varargin{1};           % case-study number
    auxdata.Pen = varargin{2};     % penalty value
    GF = varargin{3};              % geometric (size) factor for std
    auxdata.idx = varargin{4};     % size factor for idx
end


% General problem data with index -1
[opts,auxdata]=SimpleSASA_UCCD_Opts(-1,auxdata,opts);

% Case-specific problem data
auxdata.case = pcase;
[opts,auxdata]=SimpleSASA_UCCD_Opts(pcase,auxdata,opts);


% auxdata.CtrlWindow = CtrlWindow;
% auxdata.CtrlHorizon   = floor((opts.dt.nt)/(auxdata.CtrlWindow-1));

% Generate samples for MCS or gPC
switch upper(opts.UP.method)
    case 'MCS'
        % MCS-specific settings
        auxdata = SimpleSASA_smpl(auxdata,opts);
    case 'GPC'
        % gPC-specific settigns
        [auxdata, opts] = gpcopts(auxdata, opts);
end


% Solve the problem
nt = opts.dt.nt;
tic
sol = OuterLoop(auxdata,opts);
sol.time_total = toc;


% Save results with special names
if pcase == 4
    savename = strcat("sol","_",opts.UP.ctrl,"_",opts.UP.arch,"_",...
        opts.UP.form,"_",opts.UP.method,"Penalty",string(auxdata.Pen),".mat");
    % elseif pcase == 5
    %     savename = strcat("sol","_",opts.UP.ctrl,"_",opts.UP.arch,"_",...
    %         opts.UP.form,"_",opts.UP.method,"_MC_poly",".mat");

elseif pcase == 9

    savename = strcat("sol","_",opts.UP.ctrl,"_nt_",string(nt));
else
    savename = strcat("sol","_",opts.UP.ctrl,"_",opts.UP.arch,"_",...
        opts.UP.form,"_",opts.UP.method);
end


% Create folder names
Foldername = strcat('Case', string(pcase));

% save solution
if pcase ~=6 &&  pcase ~=7  &&  pcase ~=5   
    pathsol = msavename(mfilename('fullpath'),strcat('Solutions/',Foldername)); 
    save(fullfile(pathsol, savename), 'sol')
end

end



% ------------------------------------------------------------------------%
% ------------------------ gPC Opoions Function --------------------------%
% ------------------------------------------------------------------------%

function [auxdata,opts] = gpcopts(auxdata, opts)


% Obtain the collocation nodes, quadrature weights, and normalization factor 
QW = cell(1,auxdata.n.UD);
for i = 1:auxdata.n.UD
    [int,a_s,wa_s] = Gaussvar(auxdata,opts,i); % Obtain normalization factors and quadrature weights
    %name_int = strcat('int',string(i));
    %name_a_s = strcat('a_s',string(i));
    %name_wa_s = strcat('wa_s',string(i));
    name_int = 'int';
    name_a_s = 'a_s';
    name_wa_s = 'wa_s';
    S = struct(name_int, int, name_a_s, a_s, name_wa_s, wa_s);
    QW{i} = S;
end




% Create the tensor product of collocation quadrature weights
[m_a,n_a,v_a] = ndgrid(QW{1,1}.a_s,QW{1,2}.a_s,QW{1,3}.a_s);
qq = [m_a(:),n_a(:),v_a(:)];
[m_w,n_w,v_w] = ndgrid(QW{1,1}.wa_s,QW{1,2}.wa_s,QW{1,3}.wa_s);
wtc = [m_w(:),n_w(:),v_w(:)];
wcprod = kron(QW{1,1}.wa_s, kron(QW{1,2}.wa_s,QW{1,3}.wa_s));
auxdata.qq = qq;
auxdata.wtc = wtc;
auxdata.wcprod = wcprod;



for i = 1:auxdata.n.UD

    if upper(opts.q_mean(i)) ~= "UNKNOWN"
        
        for kk = 0:opts.d_i(i)

            h_(:,kk+1)=(1/((2)^(kk/2)))*hermiteH(kk,(auxdata.qq(:,i)-str2double(opts.q_mean(i)))/(...
                sqrt(2*(opts.q_std(1,i))^2)))/sqrt(QW{1,i}.int(kk +1));
        end

        QW{1,i}.h = h_;

    elseif upper(opts.q_mean(i)) == "UNKNOWN"

        QW{1,i}.h = NaN;

    end

    clear h_

end

auxdata.QW = QW;

end




% ------------------------------------------------------------------------%
% ------------------------ Outerloop Function ----------------------------%
% ------------------------------------------------------------------------%

function sol = OuterLoop(auxdata, opts)

switch upper(opts.UP.arch)
    case 'N'

        lbp = auxdata.idx*auxdata.SK;
        ubp = Inf;
        rng(8920543);
        X = rand(auxdata.n.p,1);

        options = optimoptions('fmincon','Display','Iter','Algorithm','interior-point',...
            'UseParallel',true,'MaxIter',100,'OptimalityTolerance',opts.solver.tolerance,...
            'FiniteDifferenceStepSize',sqrt(eps),'FiniteDifferenceType','central');

        if auxdata.case == 0 || auxdata.case == 1 || auxdata.case == 2
            
            tic
            [X,F] = fmincon(@(x) Objective(x,auxdata,opts,[],[]),X,[],[],[],[],lbp,ubp,[],options);
            Ctime = toc;
            sol.xpopt = X;
            sol.fval = F;

        elseif auxdata.case ==5

            opts.UP.method = 'MC_WCPOLYTOPE';
            tic
            
            for i = 1:length(auxdata.Vert)
                auxdata.i = i;
                rng(8920543);
                X = rand(auxdata.n.p,1);
                X = 4.95;
                [X,F] = fmincon(@(x) Objective(x,auxdata,opts,[],[]),X,[],[],[],[],lbp,ubp,[],options);
                Ctime = toc;
                sol.xpopt = X;
                sol.fval = F;
                sol = Createsol(sol,auxdata,opts);

                Foldername = strcat('Case', string(auxdata.case));
                savename = strcat("sol_MC_N_Det_WCpoly",string(auxdata.i),".mat");
                pathsol = msavename(mfilename('fullpath'),strcat('Solutions/',Foldername)); %#ok<NASGU>
                save(fullfile(pathsol, savename), 'sol')
                Obj(i) = F;
            end

            time = toc;

        elseif auxdata.case  == 6 || auxdata.case  == 7
            
            opts.UP.method = 'MC_WCPOLYTOPE_Magnitude';
            rng(8920543);
            X = rand(auxdata.n.p,1);
            [X,F] = fmincon(@(x) Objective(x,auxdata,opts,[],[]),X,[],[],[],[],lbp,ubp,[],options);
            Ctime = toc;
            sol.xpopt = X;
            sol.fval = F;
            sol = Createsol(sol,auxdata,opts);

            Foldername = strcat('Case', string(auxdata.case));
            if auxdata.case ==6 
                s = auxdata.GF;
            elseif auxdata.case ==7
                s = auxdata.idx;
            end
            savename = strcat("sol_MC_N_Det_WCpoly",string(s),".mat");
            pathsol = msavename(mfilename('fullpath'),strcat('Solutions/',Foldername)); %#ok<NASGU>
            save(fullfile(pathsol, savename), 'sol')

        elseif auxdata.case == 9

            tic
            lbp(2) = -auxdata.umax;
            ubp(2) = auxdata.umax;
            [X,F] = fmincon(@(x) Objective(x,auxdata,opts,[],[]),X,[],[],[],[],lbp,ubp,[],options);
            Ctime = toc;
            sol.xpopt = X;
            sol.fval = F;
            sol = Createsol(sol,auxdata,opts);

        end

    case 'SH'

         
       if auxdata.case == 6 % Use Polytopic uncertainties with various GF
            opts.UP.method = 'MC_WCPOLYTOPE_Magnitude';
            tic
            [F,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[]);
            Ctime = toc;
            sol.F = F;
            sol.T = T;
            sol.U = U;
            sol.Y = Y;
            sol.P = P;
            sol.in = in;
            sol.opts = opts;
            sol.time = Ctime;
            Foldername = strcat('Case', string(auxdata.case));
            savename = strcat("sol_MC_SH_Det_WCpoly",string(auxdata.GF),".mat");
            pathsol = msavename(mfilename('fullpath'),strcat('Solutions/',Foldername)); %#ok<NASGU>
            save(fullfile(pathsol, savename), 'sol')
            %out = [];

        elseif auxdata.case == 7 % Use Polytopic uncertainties with various idx
            opts.UP.method = 'MC_WCPOLYTOPE_Magnitude';
            tic
            [F,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[]);
            Ctime = toc;
            sol.F = F;
            sol.T = T;
            sol.U = U;
            sol.Y = Y;
            sol.P = P;
            sol.in = in;
            sol.opts = opts;
            sol.time = Ctime;
            Foldername = strcat('Case', string(auxdata.case));
            savename = strcat("sol_MC_SH_Det_WCpoly",string(auxdata.idx),".mat");
            pathsol = msavename(mfilename('fullpath'),strcat('Solutions/',Foldername)); %#ok<NASGU>
            save(fullfile(pathsol, savename), 'sol')
            %out = [];


        else
            tic
            [X,F,eflag,outpt]  = runobjconstr(auxdata,opts);
            Ctime = toc;
            sol.xpopt = X;
            sol.fval = F;
            sol.time = Ctime;
        end

end


% Reconstruct the solution for all cases but 6,7, and 8
if auxdata.case ~=5 && auxdata.case ~=6 && auxdata.case ~=7 && auxdata.case ~=9
    sol = Createsol(sol,auxdata,opts);
    sol.time = Ctime;
end


end




% ------------------------------------------------------------------------%
% ------------------------ Sampling Function -----------------------------%
% ------------------------------------------------------------------------%

function auxdata = SimpleSASA_smpl(auxdata,opts)


switch upper(opts.UP.method)
    case 'MCS'
        
        x = rand(1,auxdata.n.UD);
        rng(x(1));
        auxdata.smplk = randn(1,opts.UP.n_mcs);
        rng(x(2));
        auxdata.smplJ = auxdata.Mu_J + randn(1,opts.UP.n_mcs)*auxdata.SJ;
        rng(x(3));
        auxdata.smplY20 = 0 + randn(1,opts.UP.n_mcs)*auxdata.SY20;
        
    case 'GPC'
        
        % Only for Hermite for now
        [a_us,wa] = GaussHermite(opts.qa); % qa points in random dimension 1 created by finding the roots of the qath Hermite polynomial
        [b_us,wb] = GaussHermite(opts.qb); % qb points in random dimension 2 created by finding the roots of the qbth Hermite polynomial
        [c_us,wc] = GaussHermite(opts.qc); % qc points in random dimension 3 created by finding the roots of the qcth Hermite polynomial
        
        % Scale a within the CCD problem
        auxdata.a = a_us.*sqrt(2*auxdata.SK^2); % Complete dimension a's scaling within the CCD problem
        auxdata.b = b_us.*sqrt(2*auxdata.SJ^2) + auxdata.Mu_J;
        auxdata.c = c_us.*sqrt(2*auxdata.SY20^2) + auxdata.Mu_Y20;
        auxdata.wa = wa./sqrt(pi);
        auxdata.wb = wb./sqrt(pi);
        auxdata.wc = wc./sqrt(pi);
        
        c = 1;
        qq = zeros(opts.Q,auxdata.n.UD);
        wtc = zeros(opts.Q,auxdata.n.UD);
        wcprod = zeros(1,opts.Q);
        for ii = 1: opts.qa
            for jj = 1: opts.qb
                for kk = 1: opts.qc
                    qq(c,:) = [auxdata.a(ii), auxdata.b(jj), auxdata.c(kk)];
                    wtc(c,:) = [auxdata.wa(ii), auxdata.wb(jj), auxdata.wc(kk)];
                    wcprod(c) = auxdata.wa(ii)*auxdata.wb(jj)*auxdata.wc(kk);
                    c = c + 1;
                end
            end
        end
        auxdata.qq = qq;
        auxdata.wtc = wtc;
        auxdata.wcprod = wcprod;
  
end

end




% ------------------------------------------------------------------------%
% ------------------------ Problem Data & Settings -----------------------%
% ------------------------------------------------------------------------%

function [opts,auxdata] = SimpleSASA_UCCD_Opts(varargin)

if nargin == 0
    opts.general.displevel = 0;
    opts.general.plotflag = 1;
    opts.dt.defects = 'TR';
    opts.dt.quadraature = 'CTR';
    opts.dt.mesh = 'ED';
    opts.dt.nt = 100;
    opts.dt.ut = 100;                           
    opts.solver.tolerance = 1e-10;
    opts.general.saveflag = 1;
    % opts.solver.display = 'iter';
    opts.general.displevel = 0;

else
    pcase = varargin{1};
    auxdata = varargin{2};
    opts = varargin{3};

    switch pcase

        case -1                      % General problem data

            if~isfield(auxdata,'GF')
                auxdata.SK = 0.2;    % Standard Deviation for Plant Variable (k)
                auxdata.SJ = 0.15;   % Standard Deviation for Fixed Problem Parameter (J)
                auxdata.SY20 = 0.03; % Standard Deviation for Initial Boundary Condition (Y0)
                auxdata.Mu_J = 1;    % Mean Value for Fixed Problem Parameter (J)
                auxdata.Mu_Y20 = 0;  % Mean Value for initial state 2 (Y20)
                auxdata.t0 = 0;      % Initial time
                auxdata.tf = 1;      % Final time
                auxdata.umax = 1;    % Control limit
                auxdata.n.u = 1;     % Number of controls
                auxdata.n.x = 2;     % Number of states
                auxdata.n.p = 1;     % Number of plants
                auxdata.n.UD = 3;    % Numer of uncertain dimensions

            else                     % General problem data for various std sizes

                auxdata.SK = 0.2*auxdata.GF;     % Standard Deviation for Plant Variable (k)
                auxdata.SJ = 0.15*auxdata.GF;    % Standard Deviation for Fixed Problem Parameter (J)
                auxdata.SY20 = 0.03*auxdata.GF;  % Standard Deviation for Initial Boundary Condition (Y0)
                auxdata.Mu_J = 1;
                auxdata.Mu_Y20 = 0;
                auxdata.t0 = 0;
                auxdata.tf = 1;
                auxdata.umax = 1;
                auxdata.n.u = 1;
                auxdata.n.x = 2;
                auxdata.n.p = 1;
                auxdata.n.UD = 3;
                auxdata.idx = 3;

            end

        case 0  % Deterministic CCD 

            opts.UP.form   = "Det";      % Deterministic
            opts.UP.arch   = "N";        % Nested
            opts.UP.ctrl   = "SC";       % single-control
            opts.UP.method = "Det";      % Deterministic
            auxdata.idx    = 0;
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-6;

        case 1 % SE-UCCD using MCS

            opts.UP.form   = "Stc";      % Stochastic
            opts.UP.arch   = "N";        % Nested
            opts.UP.ctrl   = "MC";       % Multiple cntrols
            opts.UP.method = "MCS";      % MCS
            opts.UP.n_mcs  = 10^4;       % Monte Carlo number of sample
            auxdata.idx    = 3;          % Constraint shift index


        case 2 % SE-UCCD using gPC

            opts.UP.form   = "Stc";                                             % Stochastic
            opts.UP.arch   = "N";                                               % Nested
            opts.UP.ctrl   = "MC";                                              % Multi-control
            opts.UP.method = "gpc";                                             % gPC
            auxdata.idx    = 3;                                                 % Constraint shift index

            opts.q_C         = [10,10,10];                                         % number of collocation points in q dimensions (k), (J), and (Y20), respectively
            opts.Q           = opts.q_C(1)*opts.q_C(2)*opts.q_C(3);             % Total number of collocatio nodes
            opts.q_D         = ["gaussian" "gaussian" "gaussian"];              % Distributions
            opts.d_i         = [8,8,8];                                         % 1-D polynomial order (d_i) for each random dimension   
            opts.q_mean      = ["unknown" auxdata.Mu_J auxdata.Mu_Y20];         % Mean values - variables with unknown mean should be dealt with within the optimization
            opts.q_std       = [auxdata.SK, auxdata.SJ, auxdata.SY20];          % Standard deviations
            opts.PC          = (opts.d_i(1)+1)*(opts.d_i(2)+1)*(opts.d_i(3)+1); % N-D polynomial overall order



        case 3 % WCR-UCCD-SC

            auxdata.n.p = 2;               % Increase number of plants for epigraph form
            opts.dt.defects = 'TR';
            opts.dt.quadrature = 'CTR';
            opts.dt.mesh = 'ED';
            opts.dt.nt = 100;              % Number of points for (inner-loop) DT using DTQP
            opts.dt.ut = 100;               % Number of points for (outer-loop) DSS for controls
            opts.dt.u_interp = 'linear';   % Interpolation method for controls
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-6;
            opts.general.plotflag = 1;
            opts.general.saveflag = 0;
            opts.method.form = 'nonlinearprogram';
            opts.solver.function = 'ipfmincon';
            opts.solver.maxiters = 2500;
            opts.UP.form = "Det";          % Deterministic representation of uncertainties
            opts.UP.arch = "SH";           % Direct Single-shooting
            opts.UP.ctrl = "SC";           % Single-control
            opts.UP.method = "WC";         % WCR
            auxdata.idx = 3;               % Constraint shift index


        case 4 % WCR-UCCD-SC-Penlaty (multiobjective)

            auxdata.n.p = 2;              % Increase number of plants for epigraph form
            opts.dt.defects = 'TR';
            opts.dt.quadrature = 'CTR';
            opts.dt.mesh = 'ED';
            opts.dt.nt = 100;             % Number of points for (inner-loop) DT using DTQP
            opts.dt.ut = 100;             % Number of points for (outer-loop) DSS for controls
            opts.dt.u_interp = 'linear';  % Interpolation method for controls
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-6;
            opts.general.plotflag = 1;
            opts.general.saveflag = 0;
            opts.method.form = 'nonlinearprogram';
            opts.solver.function = 'ipfmincon';
            opts.solver.maxiters = 500;
            opts.UP.form = "Det";         % Deterministic representation of uncertainties
            opts.UP.arch = "SH";          % Direct Single Shooting
            opts.UP.ctrl = "SC";          % Single-control
            opts.UP.method = "WC";        % WCR
            auxdata.SF = 10;              % Scaling factor for multiobjective function
            auxdata.idx = 3;              % Constraint shift index


        case 5 % WCR-UCCD-MC-polytopic

            opts.dt.defects = 'TR';
            opts.dt.quadrature = 'CTR';
            opts.dt.mesh = 'ED';
            opts.dt.nt = 100;                            % Number of points for (inner-loop) DT using DTQP
            opts.dt.ut = 100;                             % Number of points for (outer-loop) DSS for controls
            opts.dt.u_interp = 'linear';                 % Interpolation method for controls
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-6;
            opts.general.plotflag = 1;
            opts.general.saveflag = 0;
            %opts.method.form = 'nonlinearprogram';
            opts.solver.maxiters = 500;
            opts.UP.form = "Det";                        % Deterministic representation of uncertainties
            opts.UP.arch = "N";                          % Nested
            opts.UP.ctrl = "MC";                         % Multiple-control
            opts.UP.method = "WC";                       % WCR
            auxdata.idx = 3;                             % Constraint shift index
            auxdata.Vert = createvertices(auxdata,opts); % Create Polytope vertices

        case 6 % WCR-UCCD-MC for various sizes of std

            opts.dt.defects = 'TR';
            opts.dt.quadrature = 'CTR';
            opts.dt.mesh = 'ED';
            opts.dt.nt = 100;             % Number of points for (inner-loop) DT using DTQP
            opts.dt.ut = 100;             % Number of points for (outer-loop) DSS for controls
            opts.dt.u_interp = 'linear';  % Interpolation method for controls
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-6;
            opts.general.plotflag = 1;
            opts.general.saveflag = 0;
            opts.solver.maxiters = 500;
            opts.UP.form = "Det";       % Deterministic representation of uncertainties
            opts.UP.arch = "N";        % Direct Single Shooting
            opts.UP.ctrl = "MC";        % Multiple-control
            opts.UP.method = "WC";      % WCR
            auxdata.idx = 3;            % Constraint shift index

        case 7 % WCR-UCCD-MC for various sizes of idx (uncertainty set)

            opts.dt.defects = 'TR';
            opts.dt.quadrature = 'CTR';
            opts.dt.mesh = 'ED';
            opts.dt.nt = 100;            % Number of points for (inner-loop) DT using DTQP
            opts.dt.ut = 100;             % Number of points for (outer-loop) DSS for controls
            opts.dt.u_interp = 'linear'; % Interpolation method for controls
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-6;
            opts.general.plotflag = 1;
            opts.general.saveflag = 0;
            opts.solver.maxiters = 500;
            opts.UP.form = "Det";       % Deterministic representation of uncertainties
            opts.UP.arch = "N";        % Direct Single Shooting
            opts.UP.ctrl = "MC";        % Multiple-control
            opts.UP.method = "WC";      % WCR

        case 9

            auxdata.n.p = 2;                  % Increase number of plants for epigraph form
            opts.dt.defects = 'ZO';
            opts.dt.quadrature = 'CEF';
            opts.dt.mesh = 'ED';    
            opts.dt.nt = 10;                 % Number of points for (each inner-loop MPC problem solve through DT using DTQP)
            opts.dt.u_interp = 'linear';      % Interpolation method for controls
            opts.general.displevel = 0;
            opts.solver.tolerance = 1e-5;
            opts.general.plotflag = 0;
            opts.general.Showplot = 0;
            opts.general.saveflag = 0;
            opts.solver.maxiters = 500;
            opts.UP.form = "DET";                        % Deterministic representation of uncertainties 
            opts.UP.arch = "N";                          % Direct Single-shooting
            opts.UP.ctrl = "SC-MPC2";                    % Single-control MPC
            opts.UP.method = "MPC2";                     
            auxdata.idx = 3;                             % Constraint shift index
            auxdata.Vert = createvertices(auxdata,opts); % Create Polytope vertices
            auxdata.Ctrl_Policy = 'WCR';
            auxdata.penalty = 10^2;


    end
end
end




% ------------------------------------------------------------------------%
% --------------------- Create Vertices for Polytope ---------------------%
% ------------------------------------------------------------------------%

function out = createvertices(auxdata,~)

kset = [-auxdata.idx*auxdata.SK, auxdata.idx*auxdata.SK];
Jset = [-auxdata.idx*auxdata.SJ, auxdata.idx*auxdata.SJ];
SY20set = [-auxdata.idx*auxdata.SY20, auxdata.idx*auxdata.SY20];
out = zeros(2.^auxdata.n.UD,auxdata.n.UD);
ll= 1;
for i=1:2
    for j=1:2
        for k=1:2
            out(ll,:) =[kset(i),Jset(j),SY20set(k)]; 
            ll = ll + 1;
        end
    end
end


end



% ------------------------------------------------------------------------%
% ---------------------------    END    ----------------------------------%
% ------------------------------------------------------------------------%



