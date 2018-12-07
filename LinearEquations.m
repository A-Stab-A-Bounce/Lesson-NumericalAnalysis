%solve linear equation by direct method
%written by Zhang, in Nov. 2018
function linearEquations(a,b)
scale=size(a,2);
root=zeros(1,scale);
verify=zeros(1,scale);
for k=1:scale
    brk=k;           %a(mork,k) is the biggest in row k
    for m=k:scale
        if a(m,k)>a(brk,k) brk=m; 
        end
    end
    for n=1:scale    %extend all rows
        if n~= brk && a(n,k)~=0
        z= a(n,k)/a(brk,k);
        a(n,:) = a(n,:) - z*a(brk,:);
        b(n) = b(n) - z*b(brk);
        end
    end
    tempa=a(brk,:);a(brk,:)=a(k,:);a(k,:)=tempa;
    tempb=b(brk);b(brk)=b(k);b(k)=tempb;
end

    for k=1:scale
        root(k)=(b(k)/a(k,k));
        fprintf('x%d=%.15f,\n',k,root(k));
    end 
    for k=1:scale
        for n=1:scale
        verify(k) = verify(k) + a(k,n) * root(k);
        end
    verify(k)=verify(k) - b(k);
    end
    disp('verification:');
    disp(verify);