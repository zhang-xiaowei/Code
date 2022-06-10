% gives two level orthogonal array L_{n=2^k}(2^{n-1})
% Zhang Xiaowei  Xidian University  2007.11.13
% n             =2^k,  k>1
% xz            change? Yes 1, standard L. No 0, L.
% ap            ='L/H'? L is Latin array, H is Hadamard matrix.
% for example, ZLatinArray(8,1,'L')
function sol=ZLatinArray(n,xz,ap)
% Is n equal to 2^k?
if n>300 | (ap~='L'& ap~='H') | mod(n,2)~=0
    sol='Error! n must be even and be less than or equal to 300, and ap=L or H.';
    return;
end
for i=1:10
    if n==2^i
        k=i;
        break;
    else
        i=i+1;
    end
end
if i==11 | i==1
    sol='n must be equal to 2^k, where k>1.';
    return;
end


% generate H
H2=[1 1;1 -1];
sol=H2;
if k>=1
    for i=1:k-1
        soll(1:2^i,       1:2^i) = H2(1,1).*sol;
        soll(1:2^i,       2^i+1:2^(i+1)) = H2(1,2).*sol;
        soll(2^i+1:2^(i+1), 1:2^i) = H2(2,1).*sol;
        soll(2^i+1:2^(i+1), 2^i+1:2^(i+1)) = H2(2,2).*sol;
        
        % standard L
        if xz==1
            xzsol=[];
            p=1;
            for j=1:2^i
                for m=1:2^(i+1)
                    if sol(j,:)==soll(m,1:2^i)
                        xzsol(p,:)=soll(m,:);
                        p=p+1;
                    end
                end
                % exchange
                if xzsol(p-1,2^i+1)==1
                    temp=xzsol(p-1,:);
                    xzsol(p-1,:)=xzsol(p-2,:);
                    xzsol(p-2,:)=temp;
                end
            end
            soll=xzsol;
        end
        
        sol=soll;
    end
end
if ap=='L'
    % -1 --> 2
    for i=1:2^k
        for j=1:2^k
            if sol(i,j)==-1
                sol(i,j)=2;
            end
        end
    end
    sol=sol(:,2:2^k);
end