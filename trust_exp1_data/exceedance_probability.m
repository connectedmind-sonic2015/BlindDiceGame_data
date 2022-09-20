 %% exceedance probability
 
honesty_alpha_r =  [30024,7.8874, 1.0037];
n = 10000;
%honesty_r = drchrnd(honesty_alpha_r,n)

%[eX] = eprob(honesty_r)
T_exc_p = Dir_exc_prob(honesty_alpha_r)


dishonesty_alpha_r =  [30000, 1.8717, 1.0366];
n = 10000;
%dishonesty_r = drchrnd(dishonesty_alpha_r,n)
%[DeX] = eprob(dishonesty_r)
UN_exc_p = Dir_exc_prob(dishonesty_alpha_r)

