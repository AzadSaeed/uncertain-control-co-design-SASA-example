function  [int,a_s,wa_s]  = Gaussvar(auxdata,opts,i)




% Calculated polynomial normalization factors from Xiu's book for Gaussian Distributions
int = zeros(1,opts.d_i(i)+1);
for ii = 0:opts.d_i(i)
    int(ii+1) = factorial(ii);
end



% Only for Hermite for now
[a_us,wa_us] = GaussHermite(opts.q_C(i));  % qa points in random dimension i created by finding the roots of the qath Hermite polynomial



% To reduce computational cost, perform some of the computations outside of
% the optimization problem. This can be done for the J and \xi_{2}(t_{0})
% dimensions, but for k, it needs to be calculated inside of the
% optimization problem. The "unknonw" mean indicates that a_s must be
% sscaled later (i.e. added to the mean)
% Scale
if upper(opts.q_mean(i)) == "UNKNOWN"
    a_s= a_us.*sqrt(2*(opts.q_std(i))^2); %
else
    a_s= a_us.*sqrt(2*(opts.q_std(i))^2) + str2double(opts.q_mean(i));
end

% scale the quadrature weights
wa_s = wa_us./sqrt(pi);


end