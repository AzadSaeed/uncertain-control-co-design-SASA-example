function [A,Bu,Bz] = UpdateLinearModel(xp,auxdata,opts,i,varargin)

switch upper(opts.UP.method)

    case 'MCS'

        A = [0, 1; -(xp+auxdata.smplk(i)*auxdata.SK)/auxdata.smplJ(i), 0 ];
        Bu =[0;1/auxdata.smplJ(i)];
        Bz = 0;

    case 'GPC'

        A = [0, 1; -(auxdata.qqfinal(i,1))/auxdata.qqfinal(i,2), 0 ];
        Bu =[0;1/auxdata.qqfinal(i,2)] ;
        Bz = [0;0;0;0];


    case 'DET'

        A = [0, 1; -(xp)/auxdata.Mu_J, 0 ];
        Bu =[0;1/auxdata.Mu_J] ;
        Bz = [0;0;0;0];


end

end