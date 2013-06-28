function [x y] = samplePDF(pdf, d1, d2N, N)

if nargin == 3
    if nargout <= 1
        x = samplingPDF_withoutOOP(pdf,d1,d2N);
    else
       error('samplePDF:univararg', 'Only one vector is returned for univariate distributions'); 
    end
else
    [x y] = samplingPDF_2D(pdf, [d1 d2N], N);
end

end


function x = samplingPDF_withoutOOP(pdf, domain, N)

%% Inverse transform sampling using Chebyshev technology without OOP
% June 2013
%N = 500;
global evals
evals = 1; 
% scl = 10000;  % manually scale for now.
% probability density function.
%sigma = 1/sqrt(20000); domain = [-10 10];
%pdf = @(x) 1/sqrt(2*pi*sigma.^2)*exp(-x.^2/2/sigma.^2);

% scale pdf onto [-1,1];
map = @(t) (domain(2)-domain(1)).*(t+1)./2 + domain(1);
f = @(x) pdf(map(x));

Y = simple_constructor(f);
out = sum_unit_interval(Y);
Y = Y./out;

% cumulative density function
cout = simple_cumsum(Y);
cdf = cout;

v = simple_chebpolyval(cdf); tol = 100*eps;
idx1 = find(v > tol, 1, 'first');
idx2 = find(v < 1-tol, 1, 'last');
k = length(v)-1;
x =  sin(pi*(-k:2:k)/(2*k)).';
dnew = map([x(idx1) x(idx2)]);

% scale pdf onto [-1,1];
map = @(t) (dnew(2)-dnew(1)).*(t+1)./2 + dnew(1);
f = @(x) pdf(map(x));

Y = simple_constructor(f);
out = sum_unit_interval(Y);
Y = Y./out;

x = map( generate_random_samples( Y, N ) );

% debugging
%cpdf = chebfun(pdf, domain);
%plot(x,0,'.k'), hold on,
%plot(cpdf)
end


function x = generate_random_samples( Y, N)
% cumulative density function
cout = simple_cumsum(Y);
cdf = cout;

% generate random samples from uniform distribution
r = rand(N,1);
c = cdf;

% bisection method
a = -ones(N,1); b = ones(N,1);
while norm(b-a,inf) > 1e-14
    vals = Clenshaw_evaluate(c,(a+b)/2);
    I1 = ((vals-r)<=-1e-14); I2 = ((vals-r)>=1e-14); I3 = ~I1 & ~I2;
    a = I1.*(a+b)/2 + I2.*a + I3.*(a+b)/2;
    b = I1.*b + I2.*(a+b)/2 + I3.*(a+b)/2;
end
x = (a+b)/2;
end


function cout = simple_cumsum(Y)
Y = Y(end:-1:1);
n = length(Y);
c = [0;0;Y];                                  % obtain Cheb coeffs {c_r}
cout = zeros(n-1,1);                          % initialize vector {C_r}
cout(1:n-1) = (c(3:end-1)-c(1:end-3))./...    % compute C_(n+1) ... C_2
    (2*(n:-1:2)');
cout(n,1) = c(end) - c(end-2)/2;              % compute C_1
v = ones(1,n); v(end-1:-2:1) = -1;
cout(n+1,1) = v*cout;                         % compute C_0
end

function Y = simple_constructor(f)
% simple Chebfun constructor.
for k = 2.^(3:18)
    x = sin(pi*(-k:2:k)/(2*k)).';
    vals = f(x);
    Y = [vals(end:-1:1,:) ; vals(2:end-1,:)]; % Laurent fold in columns.
    Y = fft(Y)/(2*length(x)-2);
    Y = real(Y);
    Y = Y(1:k,:);
    idx = find(abs(Y) > 10*log2(k)*eps, 1, 'last');
    if idx < k-3
        Y = Y(1:idx);
        Y(2:end-1,:) = 2*Y(2:end-1,:);
        break;
    end
end
if k == 2^18
    error
end
end

function vals = Clenshaw_evaluate(c,x)
global evals
evals = evals + 1; 
bk1 = zeros(size(x));
bk2 = bk1;
x = 2*x;
for k = 1:size(c,1)-1
    bk = c(k) + x.*bk1 - bk2;
    bk2 = bk1;
    bk1 = bk;
end
vals = c(end) + .5*x.*bk1 - bk2;
end

function out = sum_unit_interval(c)
% Integral on the unit interval
n = length(c);
if n == 1, out = c*2; return; end
c(2:2:end) = 0;
out = [2 0 2./(1-((2:n-1)).^2)]*c;
end

function v = simple_chebpolyval(c)
% Convert coefficients to values. 
c = c(:);       % Input should be a column vector
lc = length(c);
if lc == 1, v = c; return; end

ii = 2:lc-1;
c(ii) = 0.5*c(ii);
v = [c(end:-1:1) ; c(ii)];
if isreal(c)
    v = real(ifft(v));
elseif isreal(1i*c)
    v = 1i*real(ifft(imag(v)));
else
    v = ifft(v);
end
v = (lc-1)*[ 2*v(1) ; v(ii)+v(2*lc-ii) ; 2*v(lc) ];
v = v(end:-1:1);

end
