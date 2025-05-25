clc
clear
close all
rng("default")

%% add the spUtilities to the working directory
addpath("../MATLAB/spUtilities")

%% specify the test problem directory
datadir = "./datasets/";

%% load data from tim/cor/sto files

% problemName = "ssn";

% problemName = "lands2";
% problemName = "dim4";
problemName = "pgp2";
% problemName = "ieee30_ed";
% problemName = "stormG2";

fprintf("Loading %s Problem ...\n",problemName);
[rpm,stageProbs,model] = get_problem(datadir,problemName);
if length(rpm) == 2
    rpm = rpm{2};
end
nt = model.nt;
mt = model.mt;

stageNum = length(stageProbs);

%% mean value problem
params = [];
params.outputflag = 0;      % Show output during optimization
params.InfUnbdinfo = 1;
outmean = gurobi(model,params);
fprintf("Mean value: %f.\n",outmean.objval);

opts = [];
opts.maxIt = 400;
opts.savedir = './sdlogs';
if ~isfolder(opts.savedir)
    mkdir(opts.savedir)
end

gurobi_write(model,fullfile(opts.savedir,'meanvalue.lp'));

opts.writesolution = 0;

% samplePathFile = fullfile("/Users/zyzhang/Programs/outputResults/gsddptest",problemName,"IterSample.path");
samplePathFile = "";
if isfile(samplePathFile)
    sampleIdx = readmatrix(samplePathFile,'FileType','text');
    sampleIdx = sampleIdx(:,1:end-1)+1;
    opts.samplepath = sampleIdx;
    opts.maxIt = min([opts.maxIt,size(sampleIdx,1)]);
end

% initial iterate
% n1 = length(stageProbs{1}.obj);
% opts.x1 = outmean.x(1:n1);

opts.verbose = 0;
fprintf("\n************L-shaped (single-cut)************\n")
outL = LshapeSolver(stageProbs,rpm,opts);

fprintf("\n*****Stochastic Decomposition (single-cut)*****\n")

% compute the lower bound of the recourse function
n1 = length(stageProbs{1}.obj);
model.obj(1:n1) = 0;
outmean = gurobi(model,params);
opts.constL = outmean.objval;

outSD = sdSolver(stageProbs,rpm,opts);


figure;
Selected_IDX = [1:10:(length(outSD.lb)-1),length(outSD.lb)];
plot(Selected_IDX,outSD.lb(Selected_IDX),'LineWidth',2)
hold on
plot(outL.lb,'LineWidth',2)
hold on
yline(outL.lb(end),'--k','LineWidth',1.5)
grid on
legend('SD','Lshaped','optval','Location','best','fontsize',18)
xlabel("Iteration",'fontsize',18)
ylabel("Lower Bound",'fontsize',18)
title(sprintf("%s with #scenario = %d",replace(problemName,'_',''),rpm.S),'fontsize',18)


figure;
plot(outSD.x,'o','MarkerSize',12,'LineWidth',2)
hold on
plot(outL.x,'*','MarkerSize',12,'LineWidth',2)
grid on
legend('SD','Lshaped','Location','best','fontsize',18)
xlabel("Index",'fontsize',18)
ylabel("X value",'fontsize',18)
title(sprintf("%s with #scenario = %d",replace(problemName,'_',''),rpm.S),'fontsize',18)