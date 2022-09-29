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
        sol.U = out.xpopt(1:opts.dt.ut,1);

        % Get inner-loop solution
        auxdata.X = out.xpopt;

        [Fin,T,~,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[],sol.k);
        sol.Fin = Fin;
        P(3) = Y(1,2);
        sol.in = in;
        sol.opts = opts;

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

        [Fin,T,~,Y,P,in,opts] = InnerLoop(auxdata,opts,[],[],[],[],sol.k);
        sol.Fin = Fin;
        P(3) = Y(1,2);
        sol.in = in;
        sol.opts = opts;

        % save inner-loop solution
        sol.X =  auxdata.X;
        sol.T = T;
        sol.P = P;
        sol.Y = Y;

end
end