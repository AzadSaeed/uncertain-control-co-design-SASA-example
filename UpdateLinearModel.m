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

    case 'MC_WCPOLYTOPE'

        A = [0, 1; -(xp+auxdata.Vert(i,1))/(auxdata.Mu_J+auxdata.Vert(i,2)), 0 ];
        Bu =[0;1/(auxdata.Mu_J + auxdata.Vert(i,2))] ;
        Bz = [0;0;0;0];

    case 'MC_WCPOLYTOPE_MAGNITUDE'

        A = [0, 1; -(xp+auxdata.Vert(1,1))/(auxdata.Mu_J+auxdata.Vert(1,2)), 0 ];
        Bu =[0;1/(auxdata.Mu_J + auxdata.Vert(1,2))] ;
        Bz = [0;0;0;0];

    case "PARA_STUDY"
        
        if nargin > 4
            j = varargin{1,1};
        else
            i = 1;
            j=1;
        end
        A = [0, 1; -(xp)/(auxdata.J(j,i)), 0 ];
        Bu =[0;1/(auxdata.J(j,i))] ;
        Bz = [0;0;0;0];
end

end