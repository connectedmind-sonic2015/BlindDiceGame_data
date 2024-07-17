 %% exceedance probability

alpha_r =  [3300.5, 1.5];
% alpha_r =  [3327.4, 2.3, 6.3];
n = 100000;
r = drchrnd(alpha_r,n)

[eX] = eprob(r)

% 
% P_alpha_r =  [2721.5, 646.5, 1];
% n = 100000;
% P_r = drchrnd(P_alpha_r,n)




