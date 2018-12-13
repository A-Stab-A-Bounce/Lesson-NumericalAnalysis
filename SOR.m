%successive over relaxation method
%written by Zhang, in Dec. 2018

function SOR(a,b,x)
epsi = 1e-7;
n = size(a,2);
d = diag(diag(a));l = zeros(n);u = zeros(n);
eyem = eye(n);invd = eyem;
minpace=1000;minpacew=0;
fid=fopen('C:\Users\lenovo\Desktop\result.txt','a+');

for i = 2 : n
    l(i,[1:i-1]) = a(i,[1:i-1]);
    u(n+1-i,[n+2-i:n]) = a(n+1-i,[n+2-i:n]);
end    %get d,l,u

for j = 1 : n
    invd(j,j) = 1/d(j,j);
end     %get inverse of d
norm = @(mat) max(abs(mat));

for l=1:99
    w=l/50;
    fprintf(fid,'i=%d\r\n',w*50);
    p = eyem + w * invd * l;
    invp=eyem;
    for j = 1 : n
        if j ~= 1
        for k = 1 : j-1
            if p(j-1,k) ~= 0
            tempp1 = -p(j,k) / p(j-1,k);
            p(j,k) = p(j,k) + tempp1 * p(j-1,k);
            invp(j,:) = invp(j,:) + tempp1 * invp(j-1,:);
            end
        end
    end
    invp(j,:) = invp(j,:) / p(j,j);
    p(j,j)=1;
    end     %get inverse of p1
    s = invp * ((1-w) * eyem - w * invd * u);
    f = w * invp *invd * b;
    
    x1 = zeros(n,1);pace = 0;
    while norm(x1-x) > epsi
        x1 = s * x + f;
        temp = x1;x1 = x;x = temp;
        pace = pace + 1;
    end
    for i = 1: n
        fprintf(fid,'x%d=%.15f\r\n',i,x(i));
    end
        fprintf(fid,'%d\r\n\r\n',pace);
    if pace<minpace
        minpacew=l;minpace=pace;
    end
end
fprintf(fid,'optimum w=%f , optimum pace=%d',minpacew/50,minpace);
fclose(fid);
