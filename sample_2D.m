function [x,y] = sample_2D(pdf, domain, N)

mapx = @(t) (domain(2)-domain(1)).*(t+1)./2 + domain(1);
mapy = @(t) (domain(4)-domain(3)).*(t+1)./2 + domain(3);
f = @(x,y) pdf(mapx(x),mapy(y));
f = chebfun2(f);

% Is the pdf of rank-1?
if ( numel( pivots(f) ) == 1 ) 
    C = f.cols; 
    R = f.rows;
    y = mapy( sample(C, [-1,1], N) ); 
    x = mapx( sample(R, [-1,1], N) );
    return;
end


mypdf = f./sum2(f);
g = cumsum2(mypdf);

[A1,A2,A3] = chebpolyval2(g); A3 = A3.';
Xx = A1*A2*A3(:,end);
Xy = A1(end,:)*A2*A3;
% X = CC*diag(1./g.U)*RR(end,:).';
idx1 = find(1 - Xx < sqrt(eps),1,'first');
if isempty(idx1), idx1 = size(Xx,1); end
idx2 = find(Xx > sqrt(eps),1,'first');
if isempty(idx2), idx2 = 1; end

idy1 = find(1 - Xy < sqrt(eps),1,'first');
if isempty(idy1), idy1 = size(Xy,2); end
idy2 = find(Xy > sqrt(eps),1,'first');
if isempty(idy2), idy2 = 1; end
kx = size(Xx,1)-1;
xxx =  sin(pi*(-kx:2:kx)/(2*kx)).';
ky = size(Xy,2)-1;
yyy =  sin(pi*(-ky:2:ky)/(2*ky)).';
if xxx(idx1) < xxx(idx2) || yyy(idy1) < yyy(idy2)
    dnew = domain;
else
    dnew = [mapx([xxx(idx2) xxx(idx1)]) mapy([yyy(idy2) yyy(idy1)])];
end
mapx = @(t) (dnew(2)-dnew(1)).*(t+1)./2 + dnew(1);
mapy = @(t) (dnew(4)-dnew(3)).*(t+1)./2 + dnew(3);
f = @(x,y) pdf(mapx(x),mapy(y));
f = chebfun2(f);
mypdf = f./sum2(f);

margx = sum(mypdf, 1);               % marginal pdf of X
c = margx.coeffs; %c = c(end:-1:1);
x = generate_random_samples(c, [-1 1], N);

R = mypdf.rows; C = mypdf.cols; U = mypdf.pivotValues;
% Slice = R(:, x);
Slice = feval(R, x).';
C = C.values; 
% if min(size(Slice)) == 1 && N>1
%     Slice = Slice(:).';
% end
cheb_slices = C*diag(1./U)*Slice;
cheb_slices_coeffs = chebfft(cheb_slices);
sum_slices = sum_unit_interval(cheb_slices_coeffs);
scaled_slices = cheb_slices_coeffs*spdiags(1./sum_slices.',0,size(Slice,2),size(Slice,2));
scaled_slices = scaled_slices(end:-1:1,:);
idx = find(max(abs(scaled_slices),[],2)>10*eps, 1,'last');
scaled_slices = scaled_slices(1:idx,:);
y = generate_random_samples_vec(scaled_slices, [-1 1]);
x = mapx(x); y = mapy(y);

% For debugging:
% pdf = chebfun2(pdf,domain);
% contour(pdf,.01:.01:max2(pdf)), hold on,
% plot(x, y, '.k', 'markersize', 10)

end

function x = generate_random_samples( Y , dnew , N)
% cumulative density function
cout = simple_cumsum(Y);
cdf = (dnew(2)-dnew(1))/2.*cout;

% generate random samples from uniform distribution
r = rand(N,1);
c = cdf;

% bisection method
a = -ones(N,1); b = ones(N,1);
while norm(b-a,inf) > 1e-10
    vals = Clenshaw_evaluate(c,(a+b)/2);
    I1 = ((vals-r)<=-1e-14); I2 = ((vals-r)>=1e-14); I3 = ~I1 & ~I2;
    a = I1.*(a+b)/2 + I2.*a + I3.*(a+b)/2;
    b = I1.*b + I2.*(a+b)/2 + I3.*(a+b)/2;
end
x = (a+b)/2;
end
function x = generate_random_samples_vec( Y , dnew)
% generate one sample for many different pdfs.

% cumulative density function
cout = simple_cumsum_vec(Y);
cdf = (dnew(2)-dnew(1))/2.*cout;

% generate random samples from uniform distribution
N = size(Y, 2);
r = rand(N, 1);
c = cdf;

% bisection method
a = -ones(N, 1); b = ones(N, 1);
while norm(b-a,inf) > 1e-10
    vals = Clenshaw_evaluate_vec(c,(a+b)/2);
    I1 = ((vals-r)<=-1e-14); I2 = ((vals-r)>=1e-14); I3 = ~I1 & ~I2;
    a = I1.*(a+b)/2 + I2.*a + I3.*(a+b)/2;
    b = I1.*b + I2.*(a+b)/2 + I3.*(a+b)/2;
end
x = (a+b)/2;
end
function cout = simple_cumsum_vec(Y)
Y = Y(end:-1:1,:);
n = size(Y,1); m = size(Y,2);
c = [zeros(1,m);zeros(1,m);Y];                                  % obtain Cheb coeffs {c_r}
cout = zeros(n-1,m);                          % initialize vector {C_r}
% cout(1:n-1,:) = diag(1./(2*(n:-1:2)))*(c(3:end-1,:)-c(1:end-3,:));
cout(1:n-1,:) = spdiags(1./(2*(n:-1:2)).',0,n-1,n-1)*(c(3:end-1,:)-c(1:end-3,:));
cout(n,:) = c(end,:) - c(end-2,:)/2;              % compute C_1
v = ones(1,n); v(end-1:-2:1) = -1;
cout(n+1,:) = v*cout;                         % compute C_0
end

function cout = simple_cumsum(Y)
Y = Y(:); Y = Y(end:-1:1);
n = length(Y);
c = [0;0;Y];                                  % obtain Cheb coeffs {c_r}
cout = zeros(n-1,1);                          % initialize vector {C_r}
cout(1:n-1) = (c(3:end-1)-c(1:end-3))./...    % compute C_(n+1) ... C_2
    (2*(n:-1:2)');
cout(n,1) = c(end) - c(end-2)/2;              % compute C_1
v = ones(1,n); v(end-1:-2:1) = -1;
cout(n+1,1) = v*cout;                         % compute C_0
end


function vals = Clenshaw_evaluate(c,x)

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


function vals = Clenshaw_evaluate_vec(c,x)

bk1 = zeros(size(x));
bk2 = bk1;
x = 2*x;
for k = 1:size(c,1)-1
    bk = c(k,:).' + x.*bk1 - bk2;
    bk2 = bk1;
    bk1 = bk;
end
vals = c(end,:).' + .5*x.*bk1 - bk2;
end

function out = sum_unit_interval(c)
% Integral in the unit interval
n = size(c,1); c = flipud(c);
if n == 1, out = c*2; return; end
c(2:2:end,:) = 0;
out = [2 0 2./(1-((2:n-1)).^2)]*c;
end


function Y = chebfft(X)
% Vectorised chebfun fft. Converts values to coefficients in vectorised
% way.
n = size(X,1);
Y = [X(end:-1:1,:) ; X(2:end-1,:)]; % Laurent fold in columns.
if isreal(X)
    Y = fft(Y,[],1)/(2*n-2);
    Y = real(Y);
elseif isreal(1i*X)
    Y = fft(imag(Y),[],1)/(2*n-2);
    Y = 1i*real(Y);
else
    Y = fft(Y,[],1)/(2*n-2);
end
Y = Y(n:-1:1,:);
if (n > 2), Y(2:end-1,:) = 2*Y(2:end-1,:); end
end