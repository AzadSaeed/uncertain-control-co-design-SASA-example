function sol = Createsol(out, auxdata,opts)

switch upper(opts.UP.method)

    case 'GPC'

        xp = out.xpopt;

        % Scale the grid points for dimension k
        qa_grid = auxdata.qq(:,1)+ xp;
        auxdata.qqfinal = [qa_grid,auxdata.qq(:,2),auxdata.qq(:,3)];
        U = zeros(opts.dt.nt,opts.Q);
        Y1 = zeros(opts.dt.nt,opts.Q);
        Y2 = zeros(opts.dt.nt,opts.Q);
        O = zeros(1,opts.Q);

        for i=1:opts.Q
            [A,Bu,Bz] = UpdateLinearModel(out.xpopt,auxdata,opts,i);
            [o,T,u,y,P,in,Problemopts] = InnerLoop(auxdata,opts,A,Bu,Bz,i,out.xpopt);
            U(:,i) = u(:,1);
            Y1(:,i) = y(:,1);
            Y2(:,i) = y(:,2);
            O(i) = o;
        end

        % Build Polynomial bases: One-dimensional Orthogonal Polynomial bases
        h_ = zeros(opts.Q,opts.d_i(1)+1);
        for jj = 0:opts.d_i(1)
            h_(:,jj+1)=(1/((2)^(jj/2)))*hermiteH(jj,(auxdata.qqfinal(:,1)-xp)/(...
                sqrt(2*auxdata.SK^2)))/sqrt(auxdata.QW{1,1}.int(jj +1));
        end
        auxdata.QW{1,1}.h = h_;

        for ii=1:opts.Q
            m=1;
            for jj = 0:opts.d_i(1)
                for kk = 0:opts.d_i(2)
                    for ll = 0:opts.d_i(3)
                        if jj+kk+ll <= opts.PC
                            Phi(ii,m)=auxdata.QW{1,1}.h(ii,jj+1)*auxdata.QW{1,2}.h(ii,kk+1)*auxdata.QW{1,3}.h(ii,ll+1);
                            m = m +1;
                        end
                    end
                end
            end
        end
        zo = zeros(1,opts.Q);
        zY1 = zeros(length(T(:,1)),opts.Q);
        zY2 = zeros(length(T(:,1)),opts.Q);
        zU = zeros(length(T(:,1)),opts.Q);
        zhato = zeros(1,opts.PC);
        zhatU = zeros(opts.dt.nt,opts.PC);
        zhatY1 = zeros(opts.dt.nt,opts.PC);
        zhatY2 = zeros(opts.dt.nt,opts.PC);

        for ii=1:opts.PC
            for jj = 1:opts.Q
                tt = Phi(jj,ii)*auxdata.wcprod(jj);
                zo(jj) = O(jj).*tt;
                zU(:,jj) = U(:,jj).*tt;
                zY1(:,jj) = Y1(:,jj).*tt;
                zY2(:,jj) = Y2(:,jj).*tt;
            end
            zhato(:,ii) = sum(zo);
            zhatU(:,ii) = sum(zU,2);
            zhatY1(:,ii) = sum(zY1,2);
            zhatY2(:,ii) = sum(zY2,2);

        end
        sol.EU = zhatU(:,1);
        sol.EY1 = zhatY1(:,1);
        sol.EY2 = zhatY2(:,1);
        sol.EO = zhato(1);
        sol.T  = T;

        sol.VarY1 = sum(zhatY1(:,2:end).^2,2);
        sol.VarY2 = sum(zhatY2(:,2:end).^2,2);
        sol.VarU = sum(zhatU(:,2:end).^2,2);
        sol.Varo = sum(zhato(2:end).^2);
        sol.out.xpopt = out.xpopt;
        sol.p=P;
        sol.in = in;
        sol.opts = Problemopts;
        sol.out = out;
        sol.auxdata = auxdata;
        if upper(opts.UP.form) == 'STC'
            sol.objective = sol.EO;
        elseif upper(opts.UP.form) == 'PR'
            sol.objective = opts.UP.w*sol.EO + (1-opts.UP.w)*opts.UP.nf*sol.Varo;
        end



    case 'MCS'

        xp = out.xpopt;
        F = zeros(1,opts.UP.n_mcs);
        U = zeros(opts.dt.nt,opts.UP.n_mcs);
        Y1 = zeros(opts.dt.nt,opts.UP.n_mcs);
        Y2 = zeros(opts.dt.nt,opts.UP.n_mcs);
        for i=1:opts.UP.n_mcs
            [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i);
            [o,T,u,y,P,in,Problemopts] = InnerLoop(auxdata,opts,A,Bu,Bz,i,xp);
            F(i) = o;
            U(:,i) = u;
            Y1(:,i) = y(:,1);
            Y2(:,i)=y(:,2);
        end

        o = sum(F)/opts.UP.n_mcs;
        Varo = sum((F-o).^2)/(opts.UP.n_mcs-1);

        sol.T =T;
        sol.U =U;
        sol.Y1 =Y1;
        sol.Y2 =Y2;
        sol.p = P;
        sol.in = in;
        sol.xopt = xp;
        sol.out = out;
        sol.opts = Problemopts;
        sol.auxdata = auxdata;
        sol.Varo = Varo;
        sol.EF = o;


    case 'DET'
        xp = out.xpopt;
        [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,[]);
        [o,T,u,y,P,in,Problemopts] = InnerLoop(auxdata,opts,A,Bu,Bz,[],xp);
        sol.T =T;
        sol.U =u;
        sol.Y = y;
        sol.p = P;
        sol.o = o;
        sol.xopt = xp;
        sol.out = out;
        sol.opts = Problemopts;
        sol.auxdata = auxdata;
        sol.in = in;

    case 'WC'

        % outer_loop objective
        sol.Fout = out.xpopt(end,1);
        sol.k = myscale(out.xpopt(end-1,1),auxdata);
        %sol.k = out.xpopt(end-1,1);
        sol.U = out.xpopt(1:opts.dt.ut,1);

        % Get inner-loop solution
        auxdata.X = out.xpopt;
        %         if auxdata.case ==4
        %             opts.UP.method = "SC_WCPOLYTOPE";
        %             v = auxdata.X(end,1);
        %             k = myscale(auxdata.X(end-1,1),auxdata);
        %             auxdata.k = k;
        %             t_nt = linspace(0,auxdata.tf,opts.dt.ut)';
        %             u = @(t) interp1(t_nt,sol.U,t,opts.dt.u_interp);
        %             auxdata.u = u;
        %             for i = 1:length(auxdata.Vert)
        %                 [fin,T,U,Y] = InnerLoop(auxdata,opts,[],[],[],i);
        %                 Fin(i) = fin;
        %                 Y1end(i) = Y(end,1);
        %                 Y2end(i) = Y(end,2);
        %             end
        %             [sol.Fin sol.idx] = min(Fin);
        %             [Fin,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],sol.idx);
        %             t_nt = linspace(0,auxdata.tf,opts.dt.ut)';
        %             u = @(t) interp1(t_nt,sol.U,t,opts.dt.u_interp);
        %             sol.U = u(T);

        % Using DT
        %[out,F,T,U,Y,P,in,opts,idx] = WCRpolytope(auxdata.X,auxdata,opts);
        %P(1) = auxdata.Vert(idx,1);
        %P(2) = auxdata.Vert(idx,2);

        % Using simulation
        %[F,T,Y,idx] = WCRpolytope_sim(auxdata.X,auxdata,opts);
        % P(1) = auxdata.Vert(idx,1);
        % P(2) = auxdata.Vert(idx,2);
        % P(3) = auxdata.Vert(idx,3);

        if auxdata.case == 4
            opts.UP.method = "WCRPenalty";
            v = auxdata.X(end,1);
            k = myscale(auxdata.X(end-1,1),auxdata);
            auxdata.k = k;
            [Fin,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[]);
            sol.Fin = Fin;
            sol.opts = opts;
            sol.in = in;
            t_nt = linspace(0,auxdata.tf,opts.dt.ut)';
            u = @(t) interp1(t_nt,sol.U,t,opts.dt.u_interp);
            sol.U = u(T);
        else
            [Fin,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[],sol.k);
            sol.Fin = Fin;
            P(3) = Y(1,2);
            sol.in = in;
            sol.opts = opts;
        end
        % save inner-loop solution
        sol.X =  auxdata.X;
        sol.T = T;
        sol.P = P;
        sol.Y = Y;

    case 'WCRPENALTY'

        % outer_loop objective
        sol.Fout = out.xpopt(end,1);
        % sol.k = out.xpopt(end-1,1);
        sol.k = myscale(out.xpopt(end-1,1),auxdata);
        sol.U = out.xpopt(1:opts.dt.nt,1);

        % Get inner-loop solution
        auxdata.X = out.xpopt;
        if auxdata.case ==4

            % Using DT
            %             [out,F,T,U,Y,P,in,opts,idx] = WCRpolytope(auxdata.X,auxdata,opts);
            %             P(1) = auxdata.Vert(idx,1);
            %             P(2) = auxdata.Vert(idx,2);

            % Using simulation
            %[F,T,Y,idx] = WCRpolytope_sim(auxdata.X,auxdata,opts);
            % P(1) = auxdata.Vert(idx,1);
            % P(2) = auxdata.Vert(idx,2);
            % P(3) = auxdata.Vert(idx,3);
        else
            [Fin,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[]);
            sol.Fin = Fin;
            P(3) = Y(1,2);
            sol.in = in;
            sol.opts = opts;
        end
        % save inner-loop solution
        sol.X =  auxdata.X;
        sol.T = T;
        sol.P = P;
        sol.Y = Y;


    case 'MC_WCPOLYTOPE'

        xp = out.xpopt;
        [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,auxdata.i);
        [o,T,u,y,P,in,Problemopts] = InnerLoop(auxdata,opts,A,Bu,Bz,auxdata.i,xp);
        sol.T =T;
        sol.U =u;
        sol.Y = y;
        sol.p = P;
        sol.o = o;
        sol.xopt = xp;
        sol.out = out;
        sol.opts = Problemopts;
        sol.auxdata = auxdata;
        sol.in = in;

    case 'MC_WCPOLYTOPE_MAGNITUDE'

        xp = out.xpopt;
        Vert = [-auxdata.idx*auxdata.SK,auxdata.idx*auxdata.SJ,-auxdata.idx*auxdata.SY20];
        auxdata.Vert = Vert;
        [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,[]);
        [o,T,u,y,P,in,Problemopts] = InnerLoop(auxdata,opts,A,Bu,Bz,[],xp);
        sol.T =T;
        sol.U =u;
        sol.Y = y;
        sol.p = P;
        sol.o = o;
        sol.xopt = xp;
        sol.out = out;
        sol.opts = Problemopts;
        sol.auxdata = auxdata;
        sol.in = in;

    case "MPC2"

        nt = 50; % number of points in each delta_t  interval
        opts.general.Showplot = 1;
        xp = out.xpopt(1,1);
        U0 = out.xpopt(2,1);
        F  = out.fval;
        Nt = opts.dt.nt;

        % Find states at t1
        opts.UP.method = 'DET';
        [An,Bn,~] = UpdateLinearModel(xp,auxdata,opts,[]);
        opts.UP.method = "MPC2";
        tol_t          = 0;

        t_vec                 = linspace(auxdata.t0, auxdata.tf,opts.dt.nt);
        Delta_t0              = linspace(t_vec(1),t_vec(2)-tol_t,nt);
        y0                    = [0;0];
        C                     = [1, 0; 0, 1];
        D                     = 0;
        sys                   = ss(An, Bn, C, D);
        auxdata.u_zoh         = @(t) interp1(Delta_t0,U0*ones(size(Delta_t0)),t,'previous');
        [y_actual_, t_actual_, auxdata] = SimulateModel(sys, y0, Delta_t0, auxdata);

        y1_sol(1:nt,1) = interp1(t_actual_,y_actual_(:,1),Delta_t0,'linear');
        y2_sol(1:nt,1) = interp1(t_actual_,y_actual_(:,2),Delta_t0, 'linear');
        U_sol(1:nt,1)  = auxdata.u_zoh(t_actual_);
        t_sol(1:nt,1)  = Delta_t0;


       % Plot the simulation
        if opts.general.Showplot
            auxdata = PlotFunc(t_actual_, y_actual_, [], [], auxdata, 'Nominal');
        end


        % velocity will be known at t1
        auxdata.Vert = [auxdata.Vert(1,1:2);
            auxdata.Vert(3,1:2);
            auxdata.Vert(5,1:2);
            auxdata.Vert(7,1:2)];

        % Update the linear model once for all cases
        a  = zeros(2,2,4);
        bu = zeros(2,1,4);
        for i = 1:length(auxdata.Vert)
            opts.UP.method = 'MC_WCPOLYTOPE';
            [a_,bu_,~] = UpdateLinearModel(xp,auxdata,opts,i);
            opts.UP.method = 'MPC';
            a(:,:,i) = a_;
            bu(:,:,i) = bu_;
        end


        % Create one large model for multiple plants
        % This is preferred to solving each system sequentially due to
        % robust requiremnts of having the same initial control for each
        % system

        A  = zeros(auxdata.n.x*(length(auxdata.Vert)),auxdata.n.x*(length(auxdata.Vert)));
        Bu = zeros(auxdata.n.x*(length(auxdata.Vert)),(length(auxdata.Vert)));
        for i = 1:length(auxdata.Vert)
            id0 = (i-1)*auxdata.n.x+1;
            idf = id0 + auxdata.n.x-1;
            id  = id0:idf;
            A(id,id) = a(:,:,i);
            Bu(id,i)  = bu(:,:,i);
        end

        % Calculate the number of remaining intervals
        k = 1;
        auxdata.t0k = t_vec(k+1);
        auxdata.tfk = auxdata.tf;
        n_steps = opts.dt.nt-2;

        for i = 1: n_steps

            % Initial states for the current problem
            auxdata.y1k0       = y_actual_(end,1);
            auxdata.y2k0       = y_actual_(end,2);

            opts.dt.nt = opts.dt.nt - 1;
            [F,T,U,Y,~,in,opts] = InnerLoop(auxdata,opts,A,Bu,0,[],xp);
            if isnan(F) || any(isnan(T(:))) || any(isnan(U(:))) || any(isnan(Y(:)))
                warning(in.output.message);
                test(auxdata, opts, A, Bu, xp)
                u  = nan;
                y1 = nan;
                y2 = nan;
                f  = nan;
                tt = nan;
                F = nan;
            else
                u     = U;
                y1    = [Y(:,1),Y(:,3),Y(:,5),Y(:,7)];
                y2    = [Y(:,2),Y(:,4),Y(:,6),Y(:,8)];
                f     = F;
                tt    = T;
                yy1   = @(t)interp1(T,y1,t,'linear');
                yy2   = @(t)interp1(T,y2,t,'linear');
                auxdata.u_dtqp  = @(t) interp1(T,U,t,'previous');
                t1_   = tt(2:end);

                %                 if length(t1_)==1
                %                     t1_ = [tt(1)+eps, tt(2:end)];
                %                 end

                t1_del= linspace(tt(1),tt(2)-tol_t,nt);

                % Plot the part of the solutions that affects the control
                % decision in the next step
                if opts.general.Showplot

                    if length(t1_)~=1
                        up      = auxdata.u_dtqp(t1_);
                        auxdata = PlotFunc(t1_, yy1(t1_), yy2(t1_), up, auxdata, 'Future_Intervals');
                    end

                    up      = auxdata.u_dtqp(t1_del);
                    auxdata = PlotFunc(t1_del, yy1(t1_del), yy2(t1_del), up, auxdata, 'Models');

                end


                % Assign the new robust control
                switch upper(auxdata.Ctrl_Policy)
                    case 'WCR'
                        u_new = u(1,1);
                end


                % Simuate the actual model
                Delta_t           = linspace(T(1), T(2)-eps,nt);
                y0                = y_actual_(end,:);
                auxdata.u_zoh = @(t) interp1(Delta_t,u_new*ones(size(Delta_t)),t,'previous');
                [y_actual, t_actual, auxdata] = SimulateModel(sys,y0,Delta_t, auxdata);

                if opts.general.Showplot
                    auxdata = PlotFunc(Delta_t, y_actual, [], [], auxdata, 'Nominal');
                    delete(auxdata.h1);
                    delete(auxdata.h2);
                    delete(auxdata.h3);
                end

                idx0 = length(y1_sol)+1;
                idxf = idx0 + nt - 1;
                id   = idx0:idxf;
                y1_sol(id,1) = interp1(t_actual,y_actual(:,1),Delta_t,'linear');
                y2_sol(id,1) = interp1(t_actual,y_actual(:,2),Delta_t,'linear');
                U_sol(id,1)  = auxdata.u_zoh(t_actual);
                t_sol(id,1)  = Delta_t;


                % tol = auxdata.tf - t_actual(end);
                %k = k+1;
                auxdata.t0k = T(2);
                y_actual_ = y_actual(end,:);
                % F = -y_actual(end,1);

            end
        end

        opts.dt.nt = Nt;
        F = -y_actual(end,1) + auxdata.penalty*(y_actual(end,2).*y_actual(end,2));

        sol.F  = F;
        sol.xp = out.xpopt(1,1);
        sol.U0 = out.xpopt(2,1);
        sol.t  = t_sol;
        sol.U  = U_sol;
        sol.Y1 = y1_sol;
        sol.Y2 = y2_sol;
        sol.sol = out;

        ax1 = figure(1);
        savename1 = strcat('Solutions', filesep, 'Case9', filesep,'Plots',filesep, 'MPC_Y1_nt_',string(opts.dt.nt),'.pdf');
        exportgraphics(ax1,savename1)
        savename1 = strcat('Solutions', filesep, 'Case9', filesep,'Plots',filesep, 'MPC_Y1_nt_',string(opts.dt.nt));
        saveas(ax1,savename1,'fig')

        ax2 = figure(2);
        savename2 = strcat('Solutions', filesep, 'Case9', filesep,'Plots',filesep, 'MPC_Y2_nt_',string(opts.dt.nt),'.pdf');
        exportgraphics(ax2,savename2)
        savename2 = strcat('Solutions', filesep, 'Case9', filesep,'Plots',filesep, 'MPC_Y2_nt_',string(opts.dt.nt));
        saveas(ax2,savename2,'fig')


        ax3 = figure(3);
        savename3 = strcat('Solutions', filesep, 'Case9', filesep,'Plots',filesep, 'MPC_U_nt_',string(opts.dt.nt),'.pdf');
        exportgraphics(ax3,savename3)
        savename3 = strcat('Solutions', filesep, 'Case9', filesep,'Plots',filesep, 'MPC_U_nt_',string(opts.dt.nt));
        saveas(ax3,savename3,'fig')
end
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function dydt = simple_SASA_dyn(t,y,An,Bn,U)

u = U(t);
dydt = An*y + Bn*u;

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function plotmpc(x,y, xa, ya,n)

switch string(n)
    case '1'

        pos  = [70 500];
        xlab = '$t$';
        ylab = '$x_{1}$';
        %commonFigureSetup;
        commonFigureProperties;
        hf = figure(n);
        movegui(hf,pos)
        %hf = figure; % create a new figure and save handle
        hf.Color = 'w'; % change the figure background color
        hold on % do this once!

        hp = plot(x,y,'.-',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:));
        hold on
        hp = plot(xa,ya,'.',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.red(9,:),'MarkerEdgeColor',C.red(9,:));

    case '2'

        pos  = [650 500];
        xlab = '$t$';
        ylab = '$x_{2}$';

        %commonFigureSetup;
        commonFigureProperties;
        hf = figure(n);
        movegui(hf,pos)
        %hf = figure; % create a new figure and save handle
        hf.Color = 'w'; % change the figure background color
        hold on % do this once!

        hp = plot(x,y,'.-',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:));
        hold on
        hp = plot(xa,ya,'.',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.red(9,:),'MarkerEdgeColor',C.red(9,:));

    case '3'
        
        pos = [1250 500];
        xlab = '$t$';
        ylab = '$u$';

        %commonFigureSetup;
        commonFigureProperties;
        hf = figure(n);
        movegui(hf,pos)
        %hf = figure; % create a new figure and save handle
        hf.Color = 'w'; % change the figure background color
        hold on % do this once!

        hp = stairs(x,y,'.-',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:));
        hold on
        hp = stairs(xa,ya,'.',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.red(9,:),'MarkerEdgeColor',C.red(9,:));
end



if n ==1 || n ==2
    xlim([0, 1.1])
    ylim([0, 0.8])
elseif n ==3
    xlim([0, 1.1])
    ylim([-1.1, 1.1])
end


xlabel(xlab);
ylabel(ylab);
end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function auxdata = PlotFunc(t, y1, y2, u, auxdata, PlotType)

switch PlotType

    case 'Nominal'
        t_actual = t;
        y_actual = y1;
        plotmpc([],[],t_actual,y_actual(:,1),1);
        plotmpc([],[],t_actual,y_actual(:,2),2);
        plotmpc([],[],t_actual,auxdata.u_zoh(t_actual),3);

    case 'Models'

        t1 = t;
        yy1 = y1;
        yy2 = y2;
        commonFigureProperties;
        lw = 1.2;
        figure(1)
        hold on
        plot(t1,yy1,'-','LineWidth',lw,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:))

        figure(2)
        hold on
        plot(t1,yy2,'-','LineWidth',lw,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:))

        figure(3)
        hold on
        % tc1_    = linspace(auxdata.t0k,auxdata.t0k+auxdata.delta_t-eps,auxdata.nt_sim);
        stairs(t,u,'-','LineWidth',lw,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:))

    case 'Future_Intervals'

        t2 = t;
        yy1 = y1;
        yy2 = y2;
        commonFigureProperties;
        figure(1)
        hold on
        auxdata.h1 = plot(t2,yy1,'-','LineWidth',linewidth,'markersize',12,...
            'Color',C.grey(4,:),'MarkerEdgeColor',C.grey(4,:));

        figure(2)
        hold on
        auxdata.h2 = plot(t2,yy2,'o-','LineWidth',linewidth,'markersize',3,...
            'Color',C.grey(4,:),'MarkerEdgeColor',C.grey(4,:));


        figure(3)
        hold on
        auxdata.h3 = plot(t2,u,'o-','LineWidth',linewidth,'markersize',3,...
            'Color',C.grey(4,:),'MarkerEdgeColor',C.grey(4,:));

end

end


function [y_actual,t_actual, auxdata] = SimulateModel(sys, y0, Delta_t, auxdata)

wid                   = 'Control:analysis:LsimStartTime';
warning('off',wid)
U = auxdata.u_zoh(Delta_t);
t = Delta_t;
[y_actual,t_actual] = lsim(sys,U,t,y0,'zoh');

end