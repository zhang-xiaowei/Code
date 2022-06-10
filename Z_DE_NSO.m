function sol=Z_DE_NSO(LB,UB,Method)
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
A = [];  % discard pop
T1 = 32;
T2 = 16;

Iter=1;
tic
% Initialization
POP=unifrnd(LB,UB,Popsize,dim+1);
% dim+1: fvalue
for i=1:Popsize
    POP(i,dim+1)=ZMinObjFun(POP(i,1:dim));
end
% NextPOP
NextPOP = POP;
% Best
POP_Best = POP;
% XBest
POP_XBest = POP;
% XSecondBest
POP_XSecondBest = POP;
POP_XSecondBest(dim+1) = 10000;
% SC
POP_SC = zeros(Popsize,1);
% flag
POP_flag = zeros(Popsize,1);

% global min
[~,INDEX]=min(POP(:,dim+1));
GXbest=POP(INDEX,:);



while Iter<=MIter
    POP = NextPOP;
    % choose method
switch lower(Method)
    case 'de'
        F=0.5; 
        Cr=0.9;
   case 'denso'
        F=0.7;
        Cr=0.5;
    case 'mde'
        F = unifrnd(0.3,0.6);
        Cr = unifrnd(0.2,0.8);
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

    % mutation operation
    V = zeros(Popsize,dim);
    for i=1:Popsize
        % choose rand integer
        r=randperm(Popsize,5);
        while isempty(strfind(r,i)) == 0
            r=randperm(Popsize,5);
        end
        % Mutation for mutant individual
        V(i,1:dim) = POP(r(1),1:dim)+F*(POP(r(3),1:dim)-POP(r(2),1:dim));
    end
    
    % crossover operation
    U = V;
    for i =1:Popsize
        % Crossover for trial individual
        irand=unidrnd(Popsize);
        for j=1:dim
            krand=rand();
            if krand<=Cr || j==irand
            else
                U(i,j)=POP(i,j);
            end
        end
        % Caculate fx
        U(i,dim+1)=ZMinObjFun(U(i,1:dim));
    end
    
    
    % New Selection Operator
    for i = 1:Popsize
        if strcmp(Method, 'de') || strcmp(Method, 'mde')
            if U(i,dim+1) < POP(i,1+dim)
                NextPOP(i,:) = U(i,:);
            end
        elseif strcmp(Method,'denso')
            if U(i,dim+1) < POP_Best(i,dim+1)
                POP_Best(i,1:dim) =  U(i,1:dim);
                POP_Best(i,dim+1) = U(i,dim+1);
            end
            if U(i,dim+1) < POP(i,1+dim)
                NextPOP(i,1:dim+1) = U(i,:);
                A = [A;POP(i,1:1+dim)];
                if size(A,1) > Popsize
                    A(unidrnd(Popsize),:) = [];
                end
                POP_SC(i) = 0;  % SC = 0
                POP_flag(i) = 0;  % flag = 0
            else
                
                if U(i,dim+1) < POP_XBest(i,dim+1)
                    POP_XBest(i,1:dim) =  U(i,1:dim);
                    POP_XBest(i,dim+1) = U(i,dim+1);
                end
                if U(i,dim+1) < POP_XSecondBest(i,dim+1) && U(i,dim+1) > POP_XBest(i,dim+1)
                    POP_XSecondBest(i,1:dim) =  U(i,1:dim);
                    POP_XSecondBest(i,dim+1) = U(i,dim+1);
                end
                % SC = SC+1;
                POP_SC(i) = POP_SC(i)+1;
                if POP_flag(i) == 0 && POP_SC(i) > T1
                    NextPOP(i,:) = POP_XBest(j,:);
                    POP_SC(i) = 0;
                    POP_flag(i) = 1;
                end
                if POP_flag(i) == 1 && POP_SC(i) > T2
                    NextPOP(i,:) = POP_XSecondBest(j,:);
                    POP_SC(i) = 0;
                    POP_flag(i) = 2;
                end
                if POP_flag(i) == 2 && POP_SC(i) > T2
                    NextPOP(i,:) = A(unidrnd(size(A,1)),:);
                    POP_SC(i) = 0;
                    POP_flag(i) = 3;
                end
                if POP_flag(i) == 3 && POP_SC(i) > T2
                    NextPOP(i,:) = POP_Best(i,:);
                    POP_SC(i) = 0;
                    POP_flag(i) = 0;
                end
            end
        else
            disp('nothing!')
        end
        
    end  % for
    
    
    % Update GXbest and EFV, viz. current best X and f(X)
    [GXbest1,INDEX]=min(POP(:,dim+1));
    if GXbest1 < GXbest(dim+1)
        GXbest=POP(INDEX,:);
    end
    GOV(Iter)=GXbest(1+dim);
    
    % Plot
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
