function RungeKutta(fun,solution,ori,domain)
%written by Zhang, in Dec. 2018
% Dy=fun, 'solution' is the solution(func of x) of the ode, ori=[xi,yi],
% domain=[initial final], in which initial = x1
n = 4;
z = 0 : n-1;
h = 0.1 ./ 2.^z; 
y = linspace(ori(2),ori(2),size(z,2));
err = zeros(1,n);
for i = 1 : size(z,2)
    x = domain(1) : h(i) : domain(2) - h(i);
    for j = 1 : size(x,2)
    k1 = fun(x(j) , y(i));
    k2 = fun(x(j) + 0.5 * h(i) , y(i) + 0.5 * h(i) * k1);
    k3 = fun(x(j) + 0.5 * h(i) , y(i) + 0.5 * h(i) * k2);
    k4 = fun(x(j) + h(i) , y(i) + h(i) * k3);
    y(i) = y(i) + h(i) * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
    y0 = solution(x(j));
    err(i) = abs(y(i) - y0);
    if i ~= 1
        ok(i-1) = log(err(i-1) / err(i)) / log(2);
    end
end
ok(n) = 0;
for i = 1:n
    fprintf('h=%f,  err=%f  ok=%f\n',h(i),err(i),ok(i));
end
% input:
% fun = @(x,y) - x.^2 .* y.^3;
% solution = @(x) 3/(1 + x.^3);
% ori=[0,3];
% domain=[0,1.5];