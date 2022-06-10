function sol=ZMinObjFun(x)
global TESTfun;
global evaluObj;
global dim;
% PPPP XXXX wwww for "Fletcher_Powell"
global PPPP;
global XXXX;
global wwww;
MSOL=0;
%x(1:dim)=x(1:dim)-1./[1:30];
evaluObj=evaluObj+1;
switch lower(TESTfun)
%==========================================================================
    case 'sphere'
        n=dim;
        sol=0;
        for i=1:n
            sol=sol+x(i)^2;
        end
        sol=sol+MSOL;                                                   %[-100 100]
%==========================================================================        
    case 'rosenbrock'
        n=dim;
        sol=0;
        for i=1:n-1
            sol=(100.*(x(i+1)-x(i).^2).^2+(1-x(i)).^2)+sol;  %[-5.12 5.12]
        end
        sol=sol+MSOL;
%==========================================================================
    case 'deceptive'
        a=-20;
        b=0.09;
        sol=(a/(b+(x(1)-6)^2+(x(2)-6)^2)+ x(1)^2 +x(2)^2);  
        sol=sol+MSOL;                                            %[-1100 1000]
%==========================================================================
    case 'griewank'
        n=dim;
        SUM=0;
        SUM1=1;
        for i=1:n
           SUM=SUM+(1/4000)*power(x(i),2);
        end
        for i=1:n
           SUM1=SUM1*cos(x(i)*(1/sqrt(i)));
        end
        SUM=SUM-SUM1+1;
        sol=SUM+MSOL;                                             %[-500 700]
%==========================================================================
    case 'shubert'
        n=dim;
        proSUM=1;
        for i=1:n
            SUM(i)=0;
            for j=1:5
                SUM(i)=SUM(i)+j*cos((j+1)*x(i)+j);
            end
            proSUM=proSUM*SUM(i);
        end
        sol=proSUM;
        sol=sol+MSOL;                                                   %[-10 10]
%==========================================================================
    case 'schwefel'
        n=dim;
        sol=0;
        for i=1:n;
            sol=sol-x(i)*sin(sqrt(abs(x(i))));
        end
        sol=sol+MSOL;                                                   %[-500 500]
%==========================================================================
    case 'levy'
        n=dim;
        sol=0;
        x=1+(x-1)/4;
        for i=1:n-1
            sol=(x(i)-1)^2*(1+10*(sin(pi*x(i)+1)^2))+sol;
        end
        sol=sol+sin(pi*x(1))^2+(x(n)-1)^2*(1+10*(sin(2*pi*x(n)))^2);
        sol=sol+MSOL;
                                                              %[-10 10]
%==========================================================================
    case 'levy2'
        n=dim;
        sol=0;
        y=1+(x+1)/4;
        for i=1:n-1
            sol=(y(i)-1)^2*(1+10*(sin(pi*y(i+1))^2))+sol;
        end
        sol=(pi/n)*(10*sin(pi*y(1))^2+sol+(y(n)-1)^2)+...
                                             zLevyParaU(x,10,100,4,'oobj');
        sol=sol+MSOL;
                                                              %[-50 50]
%==========================================================================
    case 'levy3'                                              %[-50 50]
        n=dim;
        sol=0;
        for i=1:n-1
            sol=sol+(x(i)-1)^2*(1+sin(3*pi*x(i+1))^2);
        end
        sol=0.1*(sin(3*pi*x(1))^2+sol+(x(n)-1)^2*(1+sin(2*pi*x(n))^2))+...
            zLevyParaU(x,5,100,4,'oobj');
        sol=sol+MSOL;
    
%==========================================================================        
    case 'rastrigin'
        n=dim;
        sol=0;
        for i=1:n
            sol=x(i)^2-10*cos(2*pi*x(i))+sol;
        end
        sol=(sol+10*n); 
        sol=sol+MSOL;                                          %[-5.12 5.12]
%==========================================================================
    case 'branin'
        sol=(x(2)-(5.1/(4*pi^2))*x(1)^2+...
            (5/pi)*x(1)-6)^2+10*(1-(1/(8*pi)))...
            *cos(x(1))+10;  
        sol=sol+MSOL;                                         %[-5 10;0 15]
%==========================================================================
    case 'sixhcb'
        sol=(4*x(1)^2-2.1*x(1)^4+(1/3)*x(1)^6+...
            x(1)*x(2)-4*x(2)^2+4*x(2)^4); 
        sol=sol+MSOL;                                               %[-5 5]
%==========================================================================
    case 'beale'
        sol=((1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+...
            x(1)*x(2)^2)^2+(2.625-x(1)...
            +x(1)*x(2)^3)^2);  
        sol=sol+MSOL;                                          %[-4.5 4.5]
%==========================================================================
    case 'gp'
        T1=x(1)+x(2)+1;
        T2=19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2;
        T3=2*x(1)-3*x(2);
        T4=18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2;
        sol=((1+T1^2*T2)*(30+T3^2*T4)); 
        sol=sol+MSOL;                                                %[-2 2]
%==========================================================================
    case 'michalewics'
        n=dim;
        sol=0;
        for i=1:n
            temp=i*power(x(i),2)/pi;
            sol=sol-sin(x(i))*power(sin(temp),20);
        end
        sol=sol+MSOL;                                                   %[0 pi]
%==========================================================================
    case 'shekel5'
        m=5;                                                %m=5 7 10
        A=[4 4 4 4;...
           1 1 1 1;...
           8 8 8 8;...
           6 6 6 6;...
           3 7 3 7;...
           2 9 2 9;...
           5 5 3 3;...
           8 1 8 1;...
           6 2 6 2;...
           7 3.6 7 3.6];
       C=[1 2 2 4 4 6 3 7 5 5]/10;
       sol=0;
       for i=1:m
           sol1=0;
           for j=1:4
               sol1=sol1+(x(j)-A(i,j))^2;
           end
           sol=sol-1/(sol1+C(i));
       end
       sol=sol+MSOL;                                                    %[0 10]
%==========================================================================
    case 'shekel7'
        m=7;                                                %m=5 7 10
        A=[4 4 4 4;...
           1 1 1 1;...
           8 8 8 8;...
           6 6 6 6;...
           3 7 3 7;...
           2 9 2 9;...
           5 5 3 3;...
           8 1 8 1;...
           6 2 6 2;...
           7 3.6 7 3.6];
       C=[1 2 2 4 4 6 3 7 5 5]/10;
       sol=0;
       for i=1:m
           sol1=0;
           for j=1:4
               sol1=sol1+(x(j)-A(i,j))^2;
           end
           sol=sol-1/(sol1+C(i));
       end
       sol=sol+MSOL;                                                    %[0 10]
%==========================================================================
    case 'shekel10'
        m=10;                                                %m=5 7 10
        A=[4 4 4 4;...
           1 1 1 1;...
           8 8 8 8;...
           6 6 6 6;...
           3 7 3 7;...
           2 9 2 9;...
           5 5 3 3;...
           8 1 8 1;...
           6 2 6 2;...
           7 3.6 7 3.6];
       C=[1 2 2 4 4 6 3 7 5 5]/10;
       sol=0;
       for i=1:m
           sol1=0;
           for j=1:4
               sol1=sol1+(x(j)-A(i,j))^2;
           end
           sol=sol-1/(sol1+C(i));
       end
       sol=sol+MSOL;                                                    %[0 10]
%==========================================================================
    case 'hartman46'
        A=[   10   3    17  3.5   1.7   8;...
            0.05  10    17  0.1     8  14;...
               3 3.5   1.7   10    17   8;...
              17   8  0.05   10   0.1  14];
        P=[ 0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;...
            0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;...
            0.2348 0.1415 0.3522 0.2883 0.3047 0.6650;...
            0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
        C=[1  1.2  3  3.2];
        sol=0;
        for i=1:4
            zexpX=0;
            for j=1:6
                zexpX=zexpX+A(i,j)*(x(j)-P(i,j))^2;
            end
            sol=sol-C(i)*exp(-zexpX);
        end
        sol=sol+MSOL;
                                                              %[0 1]
%==========================================================================
    case 'hartman43'
        A=[   3   10    30;...
            0.1   10    35;...
              3   10    30;...
            0.1   10    35];
        P=[ 0.36890 0.11700 0.26730;...
            0.46990 0.43870 0.74700;...
            0.10910 0.87320 0.55470;...
            0.03815 0.57430 0.88280];
        C=[1  1.2  3  3.2];
        sol=0;
        for i=1:4
            zexpX=0;
            for j=1:3
                zexpX=zexpX+A(i,j)*(x(j)-P(i,j))^2;
            end
            sol=sol-C(i)*exp(-zexpX);
        end
        sol=sol+MSOL;
                                                             %[0 1]
%==========================================================================
    case 'ackley'
        n=dim;                                              %[-32 32]
        sol=0;
        for i=1:n
            sol=sol+cos(2*pi*x(i));
        end
        sol1=0;
        for i=1:n
            sol1=sol1+x(i)^2;
        end
        sol=-20*exp(-0.2*sqrt(sol1/n))-exp(sol/n)+exp(1)+20;
        sol=sol+MSOL;
%==========================================================================
    case 'styblinskitang'                                   %[-5 5]
        n=dim;
        sol=0;
        for i=1:n
            sol=sol+x(i)^4-16*x(i)^2+5*x(i);
        end
        sol=sol/n;
        sol=sol+MSOL;
%==========================================================================
    case 'power'                                            %[-100 100]
        n=dim;
        sol=0;
        sol1=0;
        for i=1:n
            for j=1:i
                sol1=sol1+x(j);
            end
            sol=sol+sol1^2;
        end
        sol=sol+MSOL;
%==========================================================================
    case 'xwzhang'                                          %[-100 100]
        n=dim;
        sol=0;
        for i=1:n
            sol=sol+(x(i)-i)^2;
        end
        sol=sol+MSOL;
%==========================================================================
    case 'xwzhang1'                                         %[-100 100]
        n=dim;
        sol=0;
        for i=1:n
            sol=sol+(x(i)-i)^4;
        end
        sol=sol+MSOL;
%==========================================================================
    case 'quartic'                                          %[-1.28 1.28]
        n=dim;
        sol=0;
        for i=1:n
            sol=sol+i*x(i)^4;
        end
        sol=sol+unifrnd(0,1);
        sol=sol+MSOL;
%==========================================================================
    case 'maxjdz'
        n=dim;
        for i=1:n
            tempsol(i)=abs(x(i));
        end
        sol=max(tempsol); 
        sol=sol+MSOL;                                         % [-100 100]
%==========================================================================
    case 'sumjdz'
        n=dim;
        sol=0;
        sol1=1;
        for i=1:n
            sol=sol+abs(x(i));
            sol1=sol1*abs(x(i));
        end
        sol=sol+sol1;
        sol=sol+MSOL;                                           % [-10 10]
%==========================================================================
    case 'easom'
        sol=-cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        sol=sol+MSOL;
                                                             % [-100 100]
%==========================================================================
    case 'colville'
        sol=100*(x(2)-x(1)^2)^2+(1-x(1))^2+90*(x(4)-x(3)^2)^2+...
            (1-x(3)^2)+10.1*((x(2)-1)^2+(x(4)-1)^2)+...
            19.8*(x(2)-1)*(x(4)-1);
        sol=sol+MSOL;                                            % [-10 10]
%==========================================================================
    case 'step'
        n=dim;
        sol=0;
        for i=1:n
            sol=sol+floor(x(i)+0.5)^2;
        end
        sol=sol+MSOL;                                                  % [-100 100]
%==========================================================================
    case 'perm'
        n=dim;
        sol=0;
        for k=1:n
            sol1=0;
            for i=1:n
                sol1=sol1+(i^k+0.5)*((x(i)/i)^k-1);
            end
            sol=sol+sol1^2;
        end
        sol=sol+MSOL;                                                  % [-n n]
%==========================================================================
    case 'schaffer'                                          %[-100 100]
        n=dim;
        sol=0;
        sol=0.5+(sin(sqrt(x(1)^2+x(2)^2))^2-0.5)/(1+0.001*(x(1)^2+x(2)^2))^2;
        sol=sol+MSOL;
%==========================================================================
    case 'ex2'                                   %[0 1;0 1;1.1 1.3;0 1;0 1]
        n=dim;
        sol=0;
        f=[5 -5;3 -2;2 -1;1.5 -0.5;1.2 -0.2; 1.1 -0.1];
        for i=1:6
            sol=sol+(f(i,1)-x(1)-(x(2)/(pi*i/20)^x(3)))^2+...
                (f(i,2)-(pi*i/20)*x(4)-(x(5)/(pi*i/20)^x(3)))^2;
        end
        sol=sol+MSOL;
%==========================================================================
    case 'fletcherpowell'                             %[-pi pi]
%         sol=0;
%         if evaluObj==1
%            XXXX=floor(unifrnd(-100,100,dim,dim));
%            PPPP=floor(unifrnd(-100,100,dim,dim));
%            wwww= unifrnd(-pi,pi,1,dim);
%         end
%         for i=1:dim
%             sol1=0;
%             sol2=0;
%             for j=1:dim
%                 sol1=sol1+XXXX(i,j)*sin(wwww(j))+PPPP(i,j)*cos(wwww(j));
%             end
%             for j=1:dim
%                 sol2=sol2-XXXX(i,j)*sin(x(j))+PPPP(i,j)*cos(x(j));
%             end
%             sol=sol+(sol1-sol2)^2;
%         end

%   ¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á¡Á
     if evaluObj==1
        xxxx=repmat(x(1:dim),dim,1);
        XXXX=floor(unifrnd(-100,100,dim,dim));
        PPPP=floor(unifrnd(-100,100,dim,dim));
        wwww= repmat(unifrnd(-pi,pi,1,dim),dim,1);
        wwww=XXXX.*sin(wwww)+PPPP.*cos(wwww);
     end
     xxxx=repmat(x(1:dim),dim,1);
     sol=wwww-XXXX.*sin(xxxx)-PPPP.*cos(xxxx);
     sol=sum(sol,2).^2;
     sol=sum(sol);
%==========================================================================        
    otherwise
        disp('Invalid name of function!');
end