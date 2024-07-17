% exceedance probability

T_alpha_r =  [3000.4, 1.1, 1.5];
n = 100000;
T_r = drchrnd(T_alpha_r,n)

[eX] = eprob(T_r)


UN_alpha_r =  [2178.7, 1.1, 853.2];
n = 100000;
UN_r = drchrnd(UN_alpha_r,n)

[eX] = eprob(UN_r)
