function [X,F,eflag,outpt] = runobjconstr(varargin)

if nargin == 1 % No options supplied
    opts = [];
end
auxdata = varargin{1};
opts = varargin{2};



switch upper(opts.UP.method)

    case "WC"

        lbu = -auxdata.umax*ones(opts.dt.ut ,1);
        ubu = auxdata.umax*ones(opts.dt.ut ,1);


        lbp = [-1; -inf];
        ubp = [1; inf];


        lb = [lbu; lbp];
        ub = [ubu; ubp];

        lb = lb(:); % ensure column vectors
        ub = ub(:);

        % initialize solution
        U = load('U_IG').U;
        t_nt = linspace(0,auxdata.tf,length(U))';
        u = @(t) interp1(t_nt,U,t);
        
        %t = linspace(0,auxdata.tf,opts.dt.nt)';
        t = linspace(0,auxdata.tf,opts.dt.ut)';
        Xu = u(t);
        rng(8920543); % fix seed
        Xp = rand(auxdata.n.p,1);
        Xp(1) = -0.9857; % to make sure we have the same starting point after scaling


        X = [Xu;Xp];

        [X,F,eflag,outpt] = solveWC(X,lb,ub,auxdata,opts);

 

end

end


function [X,F,eflag,outpt] = solveWC(X,lb,ub,auxdata,opts)

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint

options = optimoptions('fmincon','Display','Iter','Algorithm','interior-point',...
    'UseParallel',true,'MaxIter',500,'OptimalityTolerance',opts.solver.tolerance,...
    'FiniteDifferenceStepSize',sqrt(eps),'FiniteDifferenceType',...
    'forward', 'MaxFunctionEvaluations',3e+9);


[X,F,eflag,outpt] = fmincon(@(x) WCObjective(x,auxdata,opts),X,[],[],[],[],lb,ub,@(x) WCConst(x,auxdata,opts),options);

    function y = WCObjective(x,auxdata,opts)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = ComputeallSASA(x,auxdata,opts);
            xLast = x;
        end
        % Now compute objective function
        y = myf;
    end

    function [c,ceq] = WCConst(x,auxdata,opts)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = ComputeallSASA(x,auxdata,opts);
            xLast = x;
        end
        % Now compute constraint function
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end


end








