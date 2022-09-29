function [F,c,ceq] = ComputeallSASA(X,auxdata,opts)

switch upper(opts.UP.arch)
    case 'N'
        xp = X(end);
    case 'SH'
        auxdata.u = X(1:opts.dt.ut,1);
        auxdata.xp = X(end,1);
        xp = auxdata.xp;
end

switch upper(opts.UP.ctrl)
    case 'SC'
        switch upper(opts.UP.method)
            case 'GPC'

                qa_grid = auxdata.qq(:,1)+xp; % Scale the grid points for dimension k
                auxdata.qqfinal = [qa_grid,auxdata.qq(:,2),auxdata.qq(:,3)];
                auxdata.h1 = zeros(opts.Q,opts.da+1);
                for jj = 0:opts.da
                    auxdata.h1(:,jj+1)=(1/((2)^(jj/2)))*hermiteH(jj,(auxdata.qqfinal(:,1)-xp)/(...
                        sqrt(2*auxdata.SK^2)))/sqrt(auxdata.int1(jj +1));
                end
                for ii=1:opts.Q
                    m=1;
                    for jj = 0:opts.da
                        for kk = 0:opts.db
                            for ll = 0:opts.dc
                                if jj+kk+ll <= opts.PC
                                    auxdata.Phi(ii,m)=auxdata.h1(ii,jj+1)*auxdata.h2(ii,kk+1)*auxdata.h3(ii,ll+1);
                                    m = m +1;
                                end
                            end
                        end
                    end
                end
                switch upper(opts.UP.arch)
                    case 'N'
                        i=[];
                        [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i);
                        o = InnerLoop(auxdata,opts,A,Bu,Bz,i);
                        F= o;
                    case 'SH'

                        % o = zeros(1,opts.Q);
                        % using DTQP
                        %parfor i=1:opts.Q
                        %    o(i) = InnerLoop(auxdata,opts,[],[],[],i);
                        %end

                        % use ODE45 for faster simulations
                        Y2end = zeros(1,opts.Q);
                        o = zeros(1,opts.Q);
                        odeopts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
                        parfor i = 1:opts.Q
                            y0 = [0 auxdata.qqfinal(i,3)];
                            k = auxdata.qqfinal(i,1);
                            J = auxdata.qqfinal(i,2);
                            [~,y] = ode45(@(t,y) SASAdyn(t,y,k,J,auxdata,opts), [auxdata.t0, auxdata.tf], y0,odeopts);
                            Y2end(1,i) = y(end,2);
                            o(i) = y(end,1);
                        end
                        zo = zeros(1,opts.Q);
                        zhato = zeros(1,opts.PC);
                        for ii=1:opts.PC
                            for jj = 1:opts.Q
                                tt = auxdata.Phi(jj,ii)*auxdata.wcprod(jj);
                                zo(jj) = o(jj).*tt;
                                if isnan(zo(jj))
                                    zo(jj) = 0;
                                end
                            end
                            zhato(:,ii) = sum(zo);
                        end

                        switch upper(opts.UP.form)
                            case 'STC'
                                Eo = zhato(1);
                                o = -Eo;
                                % calculate the obj and cons in one place
                                ceq = [];
                                c =zeros(1,opts.Q*auxdata.n.x);
                                c(1:opts.Q) = -1-(Y2end/5);
                                c(opts.Q+1:opts.Q*auxdata.n.x) = (Y2end/5)-1;

                            case 'PR'
                                Eo = -zhato(1);
                                Varo = sum(zhato(2:end).^2);
                                o = opts.UP.w*Eo + (1-opts.UP.w)*opts.UP.nf*Varo;
                                % calculate the obj and cons in one place
                                ceq = [];
                                c =zeros(1,opts.Q*auxdata.n.x);
                                c(1:opts.Q) = -1-(Y2end/5);
                                c(opts.Q+1:opts.Q*auxdata.n.x) = (Y2end/5)-1;
                        end
                        F = o;



                end

            case 'WC'
                if auxdata.case == 3
                    auxdata.X = X;
                    v = X(end,1);
                    [F,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[],xp);
                    F = -v;
                    c = v-Y(end,1);
                    ceq = [];

                elseif auxdata.case == 4
                    opts.UP.method = "WCRPENALTY";
                    auxdata.X = X;
                    v = X(end,1);
                    [F,T,U,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[]);
                    F = -v;
                    c(1) = v-Y(end,1);
                    ceq = [];

                end

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
                    o(i) = InnerLoop(auxdata,opts,A,Bu,Bz,i);
                end


                h1 = zeros(opts.Q,opts.da+1);
                for jj = 0:opts.da
                    h1(:,jj+1)=(1/((2)^(jj/2)))*hermiteH(jj,(auxdata.qqfinal(:,1)-xp)/(...
                        sqrt(2*auxdata.SK^2)))/sqrt(auxdata.int1(jj +1));
                end
                for ii=1:opts.Q
                    m=1;
                    for jj = 0:opts.da
                        for kk = 0:opts.db
                            for ll = 0:opts.dc
                                if jj+kk+ll <= opts.PC
                                    Phi(ii,m)=h1(ii,jj+1)*auxdata.h2(ii,kk+1)*auxdata.h3(ii,ll+1);
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
                switch upper(opts.UP.form)
                    case 'STC'
                        Eo = zhato(1);
                        o = Eo;
                    case 'PR'
                        Eo = zhato(1);
                        Varo = sum(zhato(2:end).^2);
                        o = opts.UP.w*Eo + (1-opts.UP.w)*opts.UP.nf*Varo;
                end
                F = o;

            case 'MCS'

                parfor i=1:opts.UP.n_mcs
                    % update matrices
                    [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i);
                    % solve inner-loop
                    o(i) = InnerLoop(auxdata,opts,A,Bu,Bz,i);
                end

                F = sum(o)/opts.UP.n_mcs;
                TF = isnan(F);
                if TF
                    o_new = o(1,~isnan(o));
                    F_new = sum(o_new)/length(o_new);
                    F = F_new;
                    opts.UP.n_mcs = length(o_new);
                end
                switch upper(opts.UP.form)
                    case 'STC'
                        o = F;
                    case 'PR'
                        Eo = F;
                        Varo = sum((o-Eo).^2)/(opts.UP.n_mcs-1);
                        o = opts.UP.w*Eo + (1-opts.UP.w)*opts.UP.nf*Varo;
                end
                F = o;

        end
end
end


