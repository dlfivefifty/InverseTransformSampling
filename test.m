% Test script: Add things to test the code. 
% June, 2013. 

Battery = {
  @(x,y) exp(-10*y.^2).*exp(-100*x.^2);
  @(x,y) exp(-100*x.^2);
  @(x,y) exp(-x.^4/2 - y.^4/2).*(x-y).^2
  @(x,y) exp(-x.^2-2*y.^2).*sech(10.*x.*y)
  @(x,y) exp(-x.^2-2*y.^2).*(x-y).^2.*sech(10.*x.*y)
};

domx_array = {
    [-1 1]
    [-1 1]
    [-7 7]
    [-5 5]
    [-3 3]
};

domy_array = {
    [-1 1]
    [-1 1]
    [-7 7]
    [-4 4]
    [-3 3]
};    

for j = 1:length(Battery)
    f = Battery{j}; 
    dom1 = domx_array{j}; 
    dom2 = domy_array{j}; 
    if ~isempty(dom2)
        [X Y] = sample(f,dom1,dom2,1000);
        g = chebfun2(f,[dom1 dom2]);
        contour(g,.01:.1:max2(g)), hold on, 
        plot(X,Y,'.','markersize',6),
        axis([dom1 dom2])
    else
        X = sample(f,dom1,100);
        g = chebfun(f,dom1); 
        plot(g), hold on, plot(X,0*X,'.')
        xlim([dom1])
    end
    pause(.3), hold off
end