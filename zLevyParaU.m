function sol=zLevyParaU(x,a,k,m,zgradindex)
global dim;
% oobj/pobj/zobj/number
% a>0  and m is evel.
xsize=size(x);
if ischar(zgradindex)==1
%==========================================================================
if zgradindex=='oobj'
    sol=0;
    for i=1:dim
        if x(i)>a
            Utemp=k*(x(i)-a)^m;
        elseif x(i)<-a
            Utemp=k*(-x(i)-a)^m;
        else
            Utemp=0;
        end
        sol=sol+Utemp;
    end
elseif zgradindex=='pobj'                 %plot
        if x>a
            Utemp=k.*(x-a).^m;
        elseif x<-a
            Utemp=k.*(-x-a).^m;
        else
            Utemp=0;
        end
        sol=Utemp;
else
    sol=[0 0];
    for i=1:dim
        if x(i,1)>a & x(i,2)>a
            Utemp=k*zpower((x(i,:)-a),m);
        elseif x(i,1)<=a & x(i,1)>=-a & x(i,2)>a
            Utemp=k*zpower((x(i,:)-a),m);
            Utemp=[min(0,Utemp(1))  max(0,Utemp(2))];
        elseif x(i,1)<=a & x(i,1)>=-a & x(i,2)<=a & x(i,2)>=-a
            Utemp=[0 0];
        elseif x(i,1)<-a & x(i,2)>a
            Utemp=k*zpower((x(i,:)-a),m);
            Utemp1=k*zpower((x(i,:)+a),m);
            Utemp=[min([0,Utemp(1),Utemp1(1)])  max([0,Utemp(2),Utemp1(2)])];
        elseif x(i,1)<-a & x(i,2)>=-a & x(i,2)<=a
            Utemp=k*zpower((x(i,:)+a),m);                  % m is evel
            Utemp=[min(0,Utemp(1))  max(0,Utemp(2))];
        elseif x(i,1)<-a & x(i,2)<-a
            Utemp=k*zpower((x(i,:)+a),m);
        end
        sol=sol+Utemp;
    end   
end
%==========================================================================
else
    sol=[0 0];
    i=zgradindex;
        if x(i,1)>a & x(i,2)>a
            Utemp=m*k*zpower((x(i,:)-a),m-1);
        elseif x(i,1)<=a & x(i,1)>=-a & x(i,2)>a
            Utemp=m*k*zpower((x(i,:)-a),m-1);
            Utemp=[min(0,Utemp(1))  max(0,Utemp(2))];
        elseif x(i,1)<=a & x(i,1)>=-a & x(i,2)<=a & x(i,2)>=-a
            Utemp=[0 0];
        elseif x(i,1)<-a & x(i,2)>a
            Utemp=m*k*zpower((x(i,:)-a),m-1);
            Utemp1=m*k*zpower((znegative(x(i,:))-a),m-1);
            Utemp=[min([0,Utemp(1),Utemp1(1)])  max([0,Utemp(2),Utemp1(2)])];
        elseif x(i,1)<-a & x(i,2)>=-a & x(i,2)<=a
            Utemp=m*k*zpower((znegative(x(i,:))-a),m-1);    % m is evel
            Utemp=[min(0,Utemp(1))  max(0,Utemp(2))];
        elseif x(i,1)<-a & x(i,2)<-a
            Utemp=m*k*zpower((znegative(x(i,:))-a),m-1);
        end
        sol=Utemp;
%==========================================================================
end