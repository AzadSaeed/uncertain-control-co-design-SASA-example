function [k] = myscale(y,auxdata)

lbp = auxdata.idx*auxdata.SK;
ubp = 20; % arbitrary large number

d = 2/(ubp-lbp);
b = -(ubp+lbp)/(ubp-lbp);

k = d\(y-b);

end
