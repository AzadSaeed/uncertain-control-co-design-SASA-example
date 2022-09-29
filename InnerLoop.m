function [F,varargout] = InnerLoop(auxdata,opts,A,Bu,Bz,i,varargin)

if ~isempty(varargin)
    xp = varargin{1,1};
end


switch upper(opts.UP.form)
    case 'STC'         
        switch upper(opts.UP.ctrl) 
            case'MC'
                if upper(opts.UP.method) == "MCS"
                    M(1).right = 5; M(1).left = 0; M(1).matrix = [-1,0]; % final state

                    idx = 1;
                    LB(idx).right = 4; % initial states
                    LBmat0 = [0,auxdata.smplY20(i)];
                    LB(idx).matrix = LBmat0;
                    UB(idx).right = 4; % initial states
                    UBmat0 = [0,auxdata.smplY20(i)];
                    UB(idx).matrix = UBmat0;

                    idx = idx +1;
                    LB(idx).right = 5; % final states
                    LBmatf = [-Inf, 0];
                    LB(idx).matrix = LBmatf;
                    UB(idx).right = 5; % final states
                    UBmatf = [auxdata.umax/xp, 0];
                    UB(idx).matrix =UBmatf;

                    idx = idx +1;
                    LB(idx).right = 1; % control
                    LB(idx).matrix = -auxdata.umax;
                    UB(idx).right = 1; % control
                    UB(idx).matrix = auxdata.umax;

                    % combine
                    setup.A = A; setup.B = Bu; setup.M = M;
                    setup.UB = UB; setup.LB = LB; setup.auxdata = auxdata;
                    setup.t0 = auxdata.t0; setup.tf = auxdata.tf;

                    % form and solve the LQDO problem
                    [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

                    % optional additional outputs
                    if nargout > 1
                        varargout = {T,U,Y,P,in,opts};
                    end

                elseif upper(opts.UP.method) == "GPC"

                    M(1).right = 5; M(1).left = 0; M(1).matrix = [-1,0]; % final state

                    idx = 1;
                    LB(idx).right = 4; % initial states
                    LBmat0 = [0,auxdata.qqfinal(i,3)];
                    LB(idx).matrix = LBmat0;
                    UB(idx).right = 4; % initial states
                    UBmat0 = [0,auxdata.qqfinal(i,3)];
                    UB(idx).matrix = UBmat0;

                    idx = idx +1;
                    LB(idx).right = 5; % final states
                    LBmatf = [-Inf, 0];
                    LB(idx).matrix = LBmatf;
                    UB(idx).right = 5; % final states
                    UBmatf = [auxdata.umax/xp, 0];
                    UB(idx).matrix =UBmatf;
                    
                    
                    idx = idx +1;
                    LB(idx).right = 1; % control
                    LB(idx).matrix = -auxdata.umax;
                    UB(idx).right = 1; % control
                    UB(idx).matrix = auxdata.umax;
                    
                    % combine
                    setup.A = A; setup.B = Bu; setup.M = M;
                    setup.UB = UB; setup.LB = LB; setup.auxdata = auxdata;
                    setup.t0 = auxdata.t0; setup.tf = auxdata.tf;
                    
                    % form and solve the LQDO problem
                    [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);
                    
                    % optional additional outputs
                    if nargout > 1
                        varargout = {T,U,Y,P,in,opts};
                    end
                end
        end
        
    case 'DET'

        switch upper(opts.UP.method)
            
            case 'DET' 

                M(1).right = 5; M(1).left = 0; M(1).matrix = [-1,0]; % final state
                
                idx = 1;
                LB(idx).right = 4; % initial states
                LBmat0 = [0,0];
                LB(idx).matrix = LBmat0;
                UB(idx).right = 4; % initial states
                UBmat0 = [0,0];
                UB(idx).matrix = UBmat0;
                
                idx = idx +1;
                LB(idx).right = 5; % final states
                LBmatf = [-Inf, 0];
                LB(idx).matrix = LBmatf;
                UB(idx).right = 5; % final states
                UBmatf = [Inf, 0];
                UB(idx).matrix = UBmatf;
                
                
                idx = idx +1;
                LB(idx).right = 1; % control
                LB(idx).matrix = -auxdata.umax;
                UB(idx).right = 1; % control
                UB(idx).matrix = auxdata.umax;
                

                % combine
                setup.A = A; setup.B = Bu; setup.M = M;
                setup.UB = UB; setup.LB = LB; setup.auxdata = auxdata;
                setup.t0 = auxdata.t0; setup.tf = auxdata.tf;
                
                % form and solve the LQDO problem
                [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

                % optional additional outputs
                if nargout > 1
                    varargout = {T,U,Y,P,in,opts};
                end

            case 'WC'

                % number of controls. states,   and parameters
                n.nu = 0; n.ny = 2; n.np = 2;

                U = auxdata.X(1:opts.dt.ut,1);
                k = myscale(auxdata.X(opts.dt.ut+1,1),auxdata);
                %v = auxdata.X(opts.dt.nt+2,1);
                t_nt = linspace(0,auxdata.tf,opts.dt.ut)';
                u = @(t) interp1(t_nt,U,t,opts.dt.u_interp);

                % system dynamics
                element.dynamics = '[y2; -(k+p1)/(J+p2)*y1 + u1/(J+p2)]';
                element.parameter_list = 'k J u1 ';
                element.parameter_values = {k auxdata.Mu_J u };

                M(1).right = 5;
                M(1).left = 0;
                M(1).matrix = [1 0];

                idx = 1;
                LB(idx).right = 3; % Design Variable constraint
                LB(idx).matrix = [-auxdata.idx*auxdata.SK, -auxdata.idx*auxdata.SJ];
                UB(idx).right = 3;
                UB(idx).matrix = [auxdata.idx*auxdata.SK, auxdata.idx*auxdata.SJ];

                idx = idx + 1;
                LB(idx).right = 4;
                LB(idx).matrix = [0,0-auxdata.idx*auxdata.SY20];
                UB(idx).right = 4;
                UB(idx).matrix = [0,0+auxdata.idx*auxdata.SY20]; % initial state


                idx = idx + 1;

                LB(idx).right = 5;
                LB(idx).matrix = [-inf -0.3];
                UB(idx).right = 5;
                UB(idx).matrix = [inf 0.3];  % final state



                % guess
                Y0 = [[0,0];[1,0]];
                P0 = [0;0];
                % U0 = [[0];[0]];
                setup.guess.X = [Y0, P0];

                % combine structures
                setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
                setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

                % form and solve the LQDO problem
                [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

                %F = -Fin-v;

                %optional additional outputs
                if nargout > 1
                    varargout = {T,U,Y,P,in,opts};
                end


            case "WCRPENALTY"

                % number of controls. states,   and parameters
                n.nu = 0; n.ny = 2; n.np = 2;

                U = auxdata.X(1:opts.dt.ut,1);
                k = myscale(auxdata.X(opts.dt.ut+1,1),auxdata);
                v = auxdata.X(end,1);
                t_nt = linspace(0,auxdata.tf,opts.dt.ut)';
                u = @(t) interp1(t_nt,U,t,opts.dt.u_interp);

                % system dynamics
                element.dynamics = '[y2; -(k+p1)/(J+p2)*y1 + u1/(J+p2)]';
                element.parameter_list = 'k J u1';
                element.parameter_values = {k auxdata.Mu_J u};

                M(1).right = 5;
                M(1).left = 0;
                M(1).matrix = (1-auxdata.Pen)*[1 0];

                M(2).right = 5;
                M(2).left = 5;
                M(2).matrix = auxdata.SF*auxdata.Pen.*[0 0; 0 1];

                idx = 1;
                LB(idx).right = 3; % Design Variable constraint
                LB(idx).matrix = [-auxdata.idx*auxdata.SK, -auxdata.idx*auxdata.SJ];
                UB(idx).right = 3;
                UB(idx).matrix = [auxdata.idx*auxdata.SK, auxdata.idx*auxdata.SJ];

                idx = idx + 1;
                LB(idx).right = 4;
                LB(idx).matrix = [0,0-auxdata.idx*auxdata.SY20];
                UB(idx).right = 4;
                UB(idx).matrix = [0,0+auxdata.idx*auxdata.SY20]; % initial state
                
                
                idx = idx + 1;
                LB(idx).right = 5;
                LB(idx).matrix = [-Inf -Inf];
                UB(idx).right = 5;
                UB(idx).matrix = [Inf Inf];  % final state

                % guess
                Y0 = [[0,0];[1,0]];
                P0 = [[0.6];[0.6]];
                % U0 = [[0];[0]];
                setup.guess.X = [Y0, P0];

                % combine structures
                setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
                setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

                % form and solve the LQDO problem
                [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

                %optional additional outputs
                if nargout > 1
                    varargout = {T,U,Y,P,in,opts};
                end


            case 'MC_WCPOLYTOPE'

                n.nu = 1; n.ny = 2; n.np = 1;

                % system dynamics
                element.dynamics = '[y2; -(p1+Deltap1)/(J+Deltap2)*y1 + u1/(J+Deltap2)]';
                element.parameter_list = 'Deltap1 J Deltap2 umax';
                element.parameter_values = {auxdata.Vert(i,1) auxdata.Mu_J auxdata.Vert(i,2) auxdata.umax};


                M(1).right = 5;
                M(1).left = 0;
                M(1).matrix = -[1 0];

                idx = 1;
                LB(idx).right = 1; 
                LB(idx).matrix = -auxdata.umax;
                UB(idx).right = 1; 
                UB(idx).matrix = auxdata.umax;

                idx = idx +1;
                LB(idx).right = 3;
                LB(idx).matrix = auxdata.idx*auxdata.SK;
                UB(idx).right = 3;
                UB(idx).matrix = Inf;

                idx = idx + 1;
                LB(idx).right = 4;
                LB(idx).matrix = [0,auxdata.Vert(i,3)];
                UB(idx).right = 4;
                UB(idx).matrix = [0,auxdata.Vert(i,3)]; % initial state

                idx = idx + 1;
                LB(idx).right = 5;
                LB(idx).matrix = [-inf 0];
                UB(idx).right = 5;
                UB(idx).matrix = [inf 0];  % final state

                % guess
                Y0 = [[0,0];[1,0]];
                P0 = [[0.6];[0.6]];
                % U0 = [[0];[0]];
                setup.guess.X = [Y0, P0];

                % combine structures
                setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
                setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

                % form and solve the LQDO problem
                [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

                %optional additional outputs
                if nargout > 1
                    varargout = {T,U,Y,P,in,opts};
                end


            case 'MC_WCPOLYTOPE_MAGNITUDE'

                Vert = [-auxdata.idx*auxdata.SK,auxdata.idx*auxdata.SJ,-auxdata.idx*auxdata.SY20];
                auxdata.Vert = Vert;

                n.nu = 1; n.ny = 2; n.np = 1;

                % system dynamics
                element.dynamics = '[y2; -(p1+Deltap1)/(J+Deltap2)*y1 + u1/(J+Deltap2)]';
                element.parameter_list = 'Deltap1 J Deltap2 umax';
                element.parameter_values = {auxdata.Vert(1,1) auxdata.Mu_J auxdata.Vert(1,2) auxdata.umax };


                M(1).right = 5;
                M(1).left = 0;
                M(1).matrix = -[1 0];

                idx = 1;
                LB(idx).right = 1; 
                LB(idx).matrix = -auxdata.umax;
                UB(idx).right = 1; 
                UB(idx).matrix = auxdata.umax;

                idx = idx +1;
                LB(idx).right = 3;
                LB(idx).matrix = auxdata.idx*auxdata.SK;
                UB(idx).right = 3;
                UB(idx).matrix = Inf;

                idx = idx + 1;
                LB(idx).right = 4;
                LB(idx).matrix = [0,auxdata.Vert(1,3)];
                UB(idx).right = 4;
                UB(idx).matrix = [0,auxdata.Vert(1,3)]; % initial state

                idx = idx + 1;
                LB(idx).right = 5;
                LB(idx).matrix = [-inf 0];
                UB(idx).right = 5;
                UB(idx).matrix = [inf 0];  % final state

                % guess
                Y0 = [[0,0];[1,0]];
                P0 = [[0.6];[0.6]];
                % U0 = [[0];[0]];
                setup.guess.X = [Y0, P0];

                % combine structures
                setup.element = element; setup.M = M; setup.UB = UB; setup.LB = LB;
                setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

                % form and solve the LQDO problem
                [T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

                %optional additional outputs
                if nargout > 1
                    varargout = {T,U,Y,P,in,opts};
                end
        end

end






end
