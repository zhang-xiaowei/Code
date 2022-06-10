function Z_DE_RUN(Method)
global TESTfun;
global dim;
global evaluObj;
global MIter;
global Popsize;
clc;
% tested function set
SetFun=struct('Fun',...
    {'sphere';'ackley';'quartic';'levy2';'levy3';'griewank';'rastrigin';'rosenbrock';'sumjdz';'maxjdz';...
        'levy';'step';'michalewics';'styblinskitang';'sixhcb';'schaffer';'beale';'gp';'shekel10';'hartman46'},...
    'Dim',{30; 30;      30;       30;     30;     30;        30;         30;          30;      30; ...
         30;    30;         30;          30;             2;        2;        2;      2;    4;         6},...
    'SD',{[-100 100];[-32 32];[-1.28 1.28];[-50 50];[-50 50];[-600 600];[-5.12 5.12];[-5 10];[-10 10];[-100 100];...
        [-10 10];[-500 500];  [0 pi];        [-5 5];       [-5 5];    [-100 100]; [-4.5 4.5];[-2 2];[0 10];[0 1]});

% parameter setting    
ZEntry=[1:10];
disp(['Function ', num2str(ZEntry) ,' are chosen.']);
RUN=10;
disp(['The number of run is set as ',num2str(RUN),'.']);
Popsize=100;
disp(['Population size is ', num2str(Popsize),'.']);
MIter=2000;
disp(['Maximal iteration is ',num2str(MIter),'.']);

% tested functions required
A=[1:size(SetFun,1)];
if ZEntry~=0
    SetFun(setdiff(A,ZEntry))=[];
end
% begin to calculate
clock
SetSol={};
A=[];
clear ZEM_Results;
A=size(SetFun,1);
for i=1:A
    TESTfun=SetFun(i,1).Fun;
    dim=SetFun(i,1).Dim;
    for j=1:RUN
        SetSol{i,1}(j,:)=Z_DE(SetFun(i,1).SD(1),SetFun(i,1).SD(2),Method);
    end
    ZEM_Results(i,1)=struct('Fun',TESTfun,'Dim',dim,'SD',SetFun(i,1).SD, 'NPop',Popsize,'MIter',MIter,'NRun',RUN,...
        'sol',SetSol{i,1});
end
save(['Results\Z_DE_',Method,'_Results'],'ZEM_Results');
clock
