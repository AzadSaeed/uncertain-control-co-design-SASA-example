function [F] = Objective(X,auxdata,opts,~,~)

switch upper(opts.UP.arch)
    case 'N'
        xp = X(end);
    case 'SH'
        auxdata.u = X(1:opts.dt.nt,1);
        auxdata.xp = X(end,1);
        xp = auxdata.xp;
end

switch upper(opts.UP.ctrl)
    case 'SC'
        switch upper(opts.UP.method)

            case 'DET' % deterministic
                [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,[]);
                F = InnerLoop(auxdata,opts,A,Bu,Bz,[],xp);
        end

    case 'MC'
        switch upper(opts.UP.method)
            case 'GPC'

                % Scale the grid points for dimension k
                qa_grid = auxdata.qq(:,1)+xp;
                auxdata.qqfinal = [qa_grid,auxdata.qq(:,2),auxdata.qq(:,3)];

                % find inner-loop solution
                o = zeros(1,opts.Q);
                parfor i=1:opts.Q
                    [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i);
                    o(i) = InnerLoop(auxdata,opts,A,Bu,Bz,i,xp);
                end


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
                zhato = zeros(1,opts.PC);

                for ii=1:opts.PC
                    for jj = 1:opts.Q
                        tt = Phi(jj,ii)*auxdata.wcprod(jj);
                        zo(jj) = o(jj).*tt;
                        if isnan(zo(jj))
                            zo(jj) = 0;
                        end
                    end
                    zhato(:,ii) = sum(zo);
                end

                Eo = zhato(1);
                o = Eo;
                F = o;

            case 'MCS'

                parfor i=1:opts.UP.n_mcs
                    [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i);
                    o(i) = InnerLoop(auxdata,opts,A,Bu,Bz,i,xp);
                end

                F = sum(o)/opts.UP.n_mcs;

            case 'MC_WCPOLYTOPE'

                i = auxdata.i;

                % Satisfy dynamics at the worse-case
                [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i);
                F = InnerLoop(auxdata,opts,A,Bu,Bz,i,xp);

                % Satisfy dynamics at the nominal values
                % pdata        = auxdata;    % dummy variable
                % pdata.Vert   = zeros(8,3); % vertex values for nominal case is 0
                %[An,Bun,Bzn] = UpdateLinearModel(xp,pdata,opts,i);
                %Fn           = InnerLoop(pdata,opts,An,Bun,Bzn,i,xp);


            case 'MC_WCPOLYTOPE_MAGNITUDE'

                Vert = [-auxdata.idx*auxdata.SK,auxdata.idx*auxdata.SJ,-auxdata.idx*auxdata.SY20];
                auxdata.Vert = Vert;
                [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,[]);
                F = InnerLoop(auxdata,opts,A,Bu,Bz,[],xp);

        end

    case 'SC-MPC2'

        % Overwrite xp (from line 5) and define U0
        xp = X(1,1);
        U0 = X(2,1);

        % Find actual states at the next interval using nominal values
        auxdata.nt_sim        = 20;
        opts.UP.method        = 'DET';
        [An,Bn,~]             = UpdateLinearModel(xp,auxdata,opts,[]);
        opts.UP.method        = "MPC2";
        tol_t                 = 0;
        t_vec                 = linspace(auxdata.t0, auxdata.tf,opts.dt.nt);
        Delta_t0              = linspace(t_vec(1),t_vec(2)-tol_t,auxdata.nt_sim);
        y0                    = [0;0];
        C                     = [1, 0; 0, 1];
        D                     = 0;
        sys                   = ss(An, Bn, C, D);
        auxdata.u_zoh         = @(t) interp1(Delta_t0,U0*ones(size(Delta_t0)),t,'previous');
        [y_actual, t_actual, auxdata] = SimulateModel(sys, y0, Delta_t0, auxdata);

        % Plot the simulation
        if opts.general.Showplot
            auxdata = PlotFunc(t_actual, y_actual, [], [], auxdata, 'Nominal');
        end


        % velocity will be known at t1; so create all possible plants using
        % for vertices of the polytope for simulation
        auxdata.Vert = [auxdata.Vert(1,1:2);
            auxdata.Vert(3,1:2);
            auxdata.Vert(5,1:2);
            auxdata.Vert(7,1:2)];


        % Update the linear model once for all cases - nominal case must be
        % included and it is the index number 5
        a  = NaN(2,2,4);
        bu = NaN(2,1,4);
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
        Nt          = opts.dt.nt;
        n_steps     = opts.dt.nt-2;

        for i = 1: n_steps

            % Initial states for the current problem
            auxdata.y1k0       = y_actual(end,1);
            auxdata.y2k0       = y_actual(end,2);

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

                t1_del= linspace(tt(1),tt(2)-tol_t,auxdata.nt_sim);

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
                Delta_t           = linspace(T(1), T(2)-eps,auxdata.nt_sim);
                y0                = y_actual(end,:);
                auxdata.u_zoh = @(t) interp1(Delta_t,u_new*ones(size(Delta_t)),t,'previous');
                [y_actual, t_actual, auxdata] = SimulateModel(sys,y0,Delta_t, auxdata);

                %                 figure(4)
                %                 hhh1 = plot(T,y2,'-ob');
                %                 hold on
                %                 YYY = yy2(t1_);
                %                 hhh2 = plot(t1_,YYY,'-*r');
                %                 hold on
                %                 hhh3 = plot(Delta_t,y_actual(:,2),'-*g');
                %                 delete(hhh1)
                %                 delete(hhh2)
                %                 delete(hhh3)

                if opts.general.Showplot
                    auxdata = PlotFunc(Delta_t, y_actual, [], [], auxdata, 'Nominal');
                    delete(auxdata.h1);
                    delete(auxdata.h2);
                    delete(auxdata.h3);
                end

                % tol = auxdata.tf - t_actual(end);
                %k = k+1;
                auxdata.t0k = T(2);
                y_actual = y_actual(end,:);
                % F = -y_actual(end,1);

            end

        end

        opts.dt.nt = Nt;
        F = -y_actual(end,1) + auxdata.penalty*(y_actual(end,2).*y_actual(end,2));

end

close all

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

        hp = plot(x,y,'-',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:));
        hold on
        hp = plot(xa,ya,'-',...
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

        hp = plot(x,y,'-',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:));
        hold on
        hp = plot(xa,ya,'-',...
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

        hp = stairs(x,y,'-',...
            'linewidth',linewidth,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:));
        hold on
        hp = stairs(xa,ya,'-',...
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

function [y_actual,t_actual, auxdata] = SimulateModel(sys, y0, Delta_t, auxdata)

wid                   = 'Control:analysis:LsimStartTime';
warning('off',wid)
U = auxdata.u_zoh(Delta_t);
t = Delta_t;
[y_actual,t_actual] = lsim(sys,U,t,y0,'zoh');

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
        % tc1_    = linspace(auxdata.t0k,auxdata.t0k+auxdata.delta_t-tol_t,auxdata.nt_sim);
        stairs(t,u,'-','LineWidth',lw,'markersize',12,...
            'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(5,:))

    case 'Future_Intervals'

        t2 = t;
        yy1 = y1;
        yy2 = y2;
        commonFigureProperties;
        figure(1)
        hold on
        auxdata.h1 = plot(t2,yy1,'o-','LineWidth',linewidth,'markersize',3,...
            'Color',C.grey(4,:),'MarkerEdgeColor',C.grey(4,:));

        figure(2)
        hold on
        auxdata.h2 = plot(t2,yy2,'o-','LineWidth',linewidth,'markersize',3,...
            'Color',C.grey(4,:),'MarkerEdgeColor',C.grey(4,:));


        figure(3)
        hold on
        auxdata.h3 = stairs(t2,u,'o-','LineWidth',linewidth,'markersize',3,...
            'Color',C.grey(4,:),'MarkerEdgeColor',C.grey(4,:));

end

end