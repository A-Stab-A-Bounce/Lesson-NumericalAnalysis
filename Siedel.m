%Gauss-Siedel method
%written by Zhang, in Nov, 2018
function Siedel(a,b,x)
epsi = 1e-7;
n = size(a,2);
d = diag(diag(a));l = zeros(n);u = zeros(n);

for i = 2 : n
    l(i,[1:i-1]) = a(i,[1:i-1]);
    u(n+1-i,[n+2-i:n]) = a(n+1-i,[n+2-i:n]);
end    %get d,l,u

dl = d + l;
invdl = eye(n);
for j = 1 : n
    if j ~= 1
    for k = 1 : j-1
        if dl(j-1,k) ~= 0
        tempdl = -dl(j,k) / dl(j-1,k);
        dl(j,k) = dl(j,k) + tempdl * dl(j-1,k);
        invdl(j,:) = invdl(j,:) + tempdl * invdl(j-1,:);
        end
    end
    end
    invdl(j,:) = invdl(j,:) / dl(j,j);
    dl(j,j)=1;
end     %get inverse of d+l

s = -invdl * u;f = invdl * b; %interation form
norm = @(mat) max(abs(mat));

x1 = zeros(n,1);pace = 0;
while norm(x1-x) > epsi 
    x1 = s * x + f;
    temp = x1;x1 = x;x = temp;
    pace = pace + 1;
end
for i = 1: n
fprintf('x%d=%.15f\n',i,x(i));
end
disp(pace);