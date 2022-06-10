function sol=Z_DE_TAGUCHI(LB,UB,Method)
% x.w.zhang@126.com
% clear all
% clc
% pack;
global TESTfun; % objective function
global evaluObj;% evaluation number
global dim;     % dimension
global Popsize; % population size
global LatinArray; % Latin array
global MIter;      % maximal iter
global EFV;        % current optimal function value
tic

% Cur off if this function runs alone,
% TESTfun='rastrigin';
% LB=-5.12;
% UB=5.12;
% dim=30;
% Popsize=30;
% MIter=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iter=1;
% DE(0) <---> Taguchi(1)
Hybrid_Weight=1;     
evaluObj=0;
% 
% 32 for 30 dim; 128 for 100 dim
if dim<=3 
    Latin_Para=4; 
elseif dim>=4 && dim<=7 
    Latin_Para=8;
elseif dim>=8 && dim<=15
    Latin_Para=16;
elseif dim>=16 && dim<=31
    Latin_Para=32;
else
    Latin_Para=128;
end
LatinArray=ZLatinArray(Latin_Para,1,'L');

% Initialization
POP=unifrnd(LB,UB,Popsize,dim+1); % POP=[x,fx]
for i=1:Popsize
    POP(i,dim+1)=ZMinObjFun(POP(i,1:dim));
end
[~,INDEX]=min(POP(:,dim+1));
% Current best individual X
GXbest=POP(INDEX,:);
% Current best objective funtion f(X)
EFV=GXbest(1+dim);
GOV=inf*ones(1,MIter);

while Iter<=MIter
    for i=1:Popsize
       
        % Numerical stability because of Q1 and Q2
        temp=max(POP(:,1+dim))-min(POP(:,1+dim));
        if temp<1e-100
            break;
        end
        
        % Yield the rand integer r1,r2 and r3.
        r=randperm(Popsize,5); 
        while isempty(strfind(r,i)) == 0
            r=randperm(Popsize,5); 
        end            
                
        Q21=(POP(r(1),1+dim)-POP(r(2),1+dim))/temp;
        F21=(POP(r(2),1:dim)-POP(r(1),1:dim)).*Q21;
                
        % choose method
        % Mutation for mutant individual
        TEMPOP=ones(1,dim+1);
        switch lower(Method)
             case 'defcr1'
                 Q31=(POP(r(1),1+dim)-POP(r(3),1+dim))/temp;
                 F31=(POP(r(3),1:dim)-POP(r(1),1:dim)).*Q31;
                 TEMPOP=POP(r(1),1:dim)+(F21+F31); 
             case 'defcr2'
                 TEMPOP=POP(r(3),1:dim)+(F21); 
             case 'defcr3'
                 Q54=(POP(r(4),1+dim)-POP(r(5),1+dim))/temp;
                 F54=(POP(r(5),1:dim)-POP(r(4),1:dim)).*Q54;
                 TEMPOP=POP(r(3),1:dim)+(F21+F54);
            otherwise
        end
        
        
        % Crossover for trial individual
        % Hybrid_Weight to balance crossover between DE and Taghchi
        if rand>Hybrid_Weight
            % DE crossover
            irand=unidrnd(Popsize);
            for j=1:dim
                krand=rand();
                if krand>=0.1 && j~=irand
                    TEMPOP(j)=POP(i,j);
                end
            end
            TEMPOP(dim+1)=ZMinObjFun(TEMPOP(1:dim));
        else
            % Taguchi method crossover
            TEMPOP=DE_Taguchi_Cross([TEMPOP;POP(i,1:dim)]);
        end
        
       % Selection
       if TEMPOP(dim+1)<POP(i,1+dim)
           POP(i,:)=TEMPOP;
       end
    end  % for
    
    % Update GXbest and EFV, viz. current best X and f(X)
    [GXbest1,INDEX]=min(POP(:,dim+1));
    if GXbest1 < GXbest(dim+1)
        GXbest=POP(INDEX,:);
        EFV=GXbest(1+dim);
    end    
    GOV(Iter)=GXbest(1+dim);
    
    % Plot
%     Results=sprintf('Iter: %d   GOV: %.12f',Iter, GOV(Iter));
%     disp(Results);
%     hold on;
%     plot(Iter,GOV(Iter),'b>','MarkerSize',5);
%     title(num2str([Iter, GOV(Iter)]));
%     xlabel(num2str(evaluObj));
%     pause(0.0001);
%     hold on;
    
    % Next
    Iter=Iter+1;
end  % while

% Output Results
EachGOV=inf*ones(1,2000);
EachGOV(1:size(GOV,2))=GOV;
sol=[evaluObj toc GXbest(1:dim) GXbest(1+dim) EachGOV];
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the taguchi method to generate a better chromosome
function sol=DE_Taguchi_Cross(pop)
% NUMBER_SPLIT is the number of experiments,where it is 2^n.
global TESTfun
global dim
global LatinArray
global EFV
Latin_Row=size(LatinArray,1);
x=inf*ones(Latin_Row,dim+2);
% caculate array and generate one offspring
for m=1:Latin_Row
    for n=1:dim
        if LatinArray(m,n)==1
            x(m,n)=pop(1,n);
        else
            x(m,n)=pop(2,n);
        end
    end
    x(m,dim+1)=ZMinObjFun(x(m,:));
    
    % Accelerate convergence
    if x(m,dim+1)~=0
        if EFV>0
            DELTA=0.3;
        else
            DELTA=1.7;
        end
        % over estimation
        %x(m,dim+2)=1/((x(m,dim+1))-DELTA*EFV)^2;
        % estimation
        x(m,dim+2)=1/((x(m,dim+1))-EFV)^2;
    else
        x(m,dim+2)=1E100;
    end
end
%
sol=inf*ones(1,dim+1);
for n=1:dim
    temp=sum(x(strfind(x(:,n)',pop(1,n)),dim+2));
    if temp>sum(x(:,dim+2))-temp
        sol(n)=pop(1,n);
    else
        sol(n)=pop(2,n);
    end
end
sol(dim+1)=ZMinObjFun(sol(1:dim));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine that whether or not temp population belongs to in the feasible
% domian
function sol=IsIn(TEMPOP,LB,UB)
global dim;
temp=ones(2,dim);
temp(1,:)=temp(1,:)*LB; % upper bound
temp(2,:)=temp(2,:)*UB; % under bound
if temp(1,:)<=TEMPOP & TEMPOP<=temp(2,:)
    sol=1; % Yes
else
    sol=0; % No
end