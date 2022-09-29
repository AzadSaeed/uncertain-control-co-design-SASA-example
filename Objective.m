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


                
        end
end
end