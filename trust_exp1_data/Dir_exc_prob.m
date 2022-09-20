
function exc_p = Dir_exc_prob(alpha)
% _
% Exceedance Probability for Dirichlet-distributed random variables
% FORMAT exc_p = Dir_exc_prob(alpha)
% alpha - a 1 x K vector with Dirichlet concentration parameters
% exc_p - a 1 x K vector with Dirichlet exceedance probabilities

% Get dimensionality
%-----------------------------------------------------------------------%
K = numel(alpha);

% Analytical computation, if bivariate Dirichlet
%-----------------------------------------------------------------------%
if K == 2
% using the Beta CDF
    exc_p(1) = 1 - betainc(1/2,alpha(1),alpha(2));
    exc_p(2) = 1 - exc_p(1);
end;

% Numerical integration, if multivariate Dirichlet
%-----------------------------------------------------------------------%
if K > 2
% using Gamma CDFs
    exc_p = zeros(1,K);
    for j = 1:K
       f = @(x) integrand(x,alpha(j),alpha([1:K]~=j));
      exc_p(j) = integral(f,0,alpha(j)) + integral(f,alpha(j),Inf);
    end;
end;
% Integrand function for numerical integration
%-----------------------------------------------------------------------%
function p = integrand(x,aj,ak)

% product of Gamma CDFs
p = ones(size(x));
for k = 1:numel(ak)
    p = p .* gammainc(x,ak(k));
end;
% times a Gamma PDF
p = p .* exp((aj-1).*log(x) - x - gammaln(aj));