%find the eigenvalue and eigenvector of a symmetrical matrix
%written by Zhang, in Dec. 2018
function eigen(a) %a , symmetrical
n = size(a,1);v = eye(n);
epsi = 1e-14;
b = a - diag(diag(a));
%format long e;
while max(b(:)) > epsi 
    m = eye(n);
    p = 1;q = 2;
    for j = 1:n
        for k = 1:n 
        if j ~= k
            if abs(b(j,k)) > abs(b(p,q))
                p = j;q = k;
            end
        end
        end
    end     %find non-diagonal num a(p,q)
    s = (a(q,q) - a(p,p)) / (2 * a(p,q));
    if s >= 0
        t = 1 / (s + sqrt(1 + s^2));
    else
        t = 1 / (s - sqrt(1 + s^2));
    end    
    c = 1 / sqrt(1 + t^2);
    d = t * c;
    m(p,p) = c;m(q,q) = c;m(p,q) = d;m(q,p) = -d;
    a = m' * a * m;
    v = v * m;
    b = a - diag(diag(a));
end
for i = 1:n
    fprintf('r%d = %.15f\n',i,a(i,i));
    fprintf('v%d = ',i);
    disp(v(:,i)');
end