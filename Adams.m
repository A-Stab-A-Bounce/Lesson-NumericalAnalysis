function Adams(fun,solution,ori,domain)
%written by Zhang, in Dec. 2018
% Dy=fun, 'solution' is the solution(func of x) of the ode, ori=[xi,yi],
% domain=[initial final], in which initial = x1
n = 4;
z = 0 : n-1;
h = 0.1 ./ 2.^z; 
yb = linspace(ori(2),ori(2),size(z,2));
yf = linspace(ori(2),ori(2),size(z,2));
err = zeros(1,n);
    
for i = 1 : size(z,2)
    x = domain(1) : h(i) : domain(2) - h(i);
    k1 = fun(x(1) , yb(i));
    k2 = fun(x(1) + 0.5 * h(i) , yb(i) + 0.5 * h(i) * k1);
    k3 = fun(x(1) + 0.5 * h(i) , yb(i) + 0.5 * h(i) * k2);
    k4 = fun(x(1) + h(i) , yb(i) + h(i) * k3);
    yf(i) = yb(i) + h(i) * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end %initialization

for i = 1 :size(z,2)   
    for j = 3 : size(x,2)
        ytemp = yf(i);
        for num = 1:100
            ytemp1 = ytemp;
            ytemp = yf(i) + h(i) * (5 * fun(x(j),ytemp) + 8 * fun(x(j-1),yf(i)) - fun(x(j-2),yb(i))) / 12;
            if abs(ytemp1 - ytemp) < 10e-5
                break;
            end
        end
        yb(i) = yf(i);yf(i) = ytemp;
    end

    y0 = solution(x(j));
    err(i) = abs(yf(i) - y0);
    if i ~= 1
        ok(i-1) = log(err(i-1) / err(i)) / log(2);
    end
end
ok(n) = 0;
for i = 1:n
    fprintf('h=%f,  err=%f  ok=%f\n',h(i),err(i),ok(i));
end
