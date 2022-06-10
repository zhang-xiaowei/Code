function sol=Z_DE(LB,UB,Method)
% clear all
% clc
% pack;
 
global TESTfun; % objective function
global evaluObj;% evaluation number
global dim;     % dimension
global Popsize; % population size
global MIter;      % maximal iter


evaluObj=0;
POP=[];  % POP=[x,fx]
GOV=inf*ones(1,MIter);  % Global Optimal Value

Iter=1;
tic
% Initialization
POP=unifrnd(LB,UB,Popsize,dim+1);
for i=1:Popsize
    POP(i,dim+1)=ZMinObjFun(POP(i,1:dim));
end
[~,INDEX]=min(POP(:,dim+1));
GXbest=POP(INDEX,:);


while Iter<=MIter
    for i=1:Popsize
        TEMPOP=(UB+1)*ones(1,dim);
              
         % choose method
         switch lower(Method)
             case 'de'
                 F=0.5;
                 Cr=0.9;
             case 'deg'
                 F=normrnd(0.5,0.3);
                 Cr=0.9;
             case 'de04'
                 F=0.4+0.4*rand;
                 Cr=0.9;
             case 'dem'
                 temp = abs(max(POP(:,1+dim))/min(POP(:,1+dim)));
                 if temp<1
                     F=max(0.4,1-temp);
                 else
                     F=max(0.4,1-1/temp);
                 end
                 Cr=0.5;
             otherwise
         end
         
         % choose rand integer
        r=randperm(Popsize,5); 
        while isempty(strfind(r,i)) == 0
            r=randperm(Popsize,5); 
        end   
        
         % Mutation for mutant individual
         TEMPOP=POP(r(1),1:dim)+F*(POP(r(3),1:dim)-POP(r(2),1:dim)); 
         
         % Crossover for trial individual
         irand=unidrnd(Popsize);
         for j=1:dim
             krand=rand();
             if krand<=Cr || j==irand
             else
                 TEMPOP(j)=POP(i,j);
             end
         end
         % Caculate fx
         TEMPOP(dim+1)=ZMinObjFun(TEMPOP(1:dim));
         % Selection
         if TEMPOP(dim+1)<POP(i,1+dim)
             POP(i,:)=TEMPOP;
         end
         
    end  % for
             

    % Update GXbest and EFV, viz. current best X and f(X)
    [GXbest1,INDEX]=min(POP(:,dim+1));
    if GXbest1 < GXbest(dim+1)
        GXbest=POP(INDEX,:);
    end    
    GOV(Iter)=GXbest(1+dim);
    
%     % Plot
    Results=sprintf('Iter: %d   GOV: %.12f',Iter, GOV(Iter));
    disp(Results);
    hold on;
    plot(Iter,GOV(Iter),'b>','MarkerSize',5);
    title(num2str([Iter, GOV(Iter)]));
    xlabel(num2str(evaluObj));
    pause(0.0001);
    hold on;

    Iter=Iter+1;
    
end  % while

% Output Results
EachGOV=inf*ones(1,2000);
EachGOV(1:size(GOV,2))=GOV;
sol=[evaluObj toc GXbest(1:dim) GXbest(1+dim) EachGOV];
