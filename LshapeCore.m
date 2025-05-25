function [stageProbs,out] = LshapeCore(stageProbs,rpm,opts,k,x1)
% Inputs:
% stageProbs{t}: the main program at stage t
% rpm: uncertainty at stage t
% k: iteration index
% opts: algorithm options
% - opts.verbose: 1 for printing the iteration information
% - opts.logdir: iterate log save path
% x1: initial iterate at stage 1
% 
% Outputs:
% out.ub: the upper bound at iteration k
% out.lb: the lower bound at iteration k
% out.x: the solution at iteration k
% out.stop: stopping flag
% 
% Author: Zhiyuan Zhang
% Date: May 25, 2025

stageNum = length(stageProbs);
params = [];
params.outputflag = 0;
params.InfUnbdinfo = 1;

x = cell(stageNum,1);

% stage 1
t = 1;

if ~isempty(x1)
    x{t} = x1;
    lb = -inf;
else
    programMain = [];
    programMain.obj = [stageProbs{t}.obj; 1.0];
    programMain.rhs = [stageProbs{t}.r; stageProbs{t}.alpha];
    programMain.A = [stageProbs{t}.D;stageProbs{t}.neg_beta_eta'];
    programMain.modelname = sprintf('It%dFPstage%dPro',k,t);
    programMain.varnames = stageProbs{t}.varnames;
    programMain.sense = strcat(stageProbs{t}.sense,repmat('>',1,length(stageProbs{t}.alpha)));

    if isfield(stageProbs{t},'ub')
        programMain.ub = stageProbs{t}.ub;
    end

    if isfield(stageProbs{t},'lb')
        programMain.lb = stageProbs{t}.lb;
        if isempty(stageProbs{t}.alpha)
            % assume eta = 0
            programMain.lb(end,1) = 0;
        end
    end
    outxt = gurobi(programMain,params);
    if strcmp(outxt.status,'OPTIMAL')
        x{t} = outxt.x(1:end-1);
        lb = outxt.objval;
    end
end

% stage 2: solve the dual to obtain the cut.
t = 2;
Cx = stageProbs{t}.C * x{t-1};
programMain = [];
programMain.A = stageProbs{t}.D;
programMain.obj = stageProbs{t}.obj;
programMain.sense = stageProbs{t}.sense;
programMain.varnames = stageProbs{t}.varnames;
if isfield(stageProbs{t},'ub') && isfield(stageProbs{t},'lb')
    programMain.lb = stageProbs{t}.lb;
    programMain.ub = stageProbs{t}.ub;
elseif isfield(stageProbs{t},'ub')
    programMain.ub = stageProbs{t}.ub;
elseif isfield(stageProbs{t},'lb')
    programMain.lb = stageProbs{t}.lb;
end


alphaS = zeros(rpm.S,1);
betaS = zeros(size(x{t-1},1), rpm.S);
objS = zeros(rpm.S,1);
for i = 1:rpm.S
    rti = stageProbs{t}.r + rpm.val(:,i);
    programMain.rhs = rti - Cx;
    outxti = gurobi(programMain,params);
    if strcmp(outxti.status,'OPTIMAL')
        objS(i,1) = outxti.objval;
        if isfield(programMain,'ub') && isfield(programMain,'lb')
            idxub = programMain.ub < inf;
            idxlb = programMain.lb > -inf;
            alphaS(i,1) = rti' * outxti.pi + programMain.ub(idxub)'*min(outxti.rc(idxub),0) + ...
                programMain.lb(idxlb)'*max(outxti.rc(idxlb),0);
        elseif isfield(programMain,'ub')
            idxub = programMain.ub < inf;
            alphaS(i,1) = rti' * outxti.pi + programMain.ub(idxub)'*min(outxti.rc(idxub),0);
        elseif isfield(programMain,'lb')
            idxlb = programMain.lb > -inf;
            alphaS(i,1) = rti' * outxti.pi + programMain.lb(idxlb)'*max(outxti.rc(idxlb),0);
        else
            alphaS(i,1) = rti' * outxti.pi;
        end
        betaS(:,i) = - stageProbs{t}.C' * outxti.pi;
    end
end

alphan = alphaS' * rpm.probs;
betan = betaS * rpm.probs;
ojbstage2 = objS' * rpm.probs;

newCutFlag = addCut2Pool([alphan;betan], [stageProbs{t-1}.alpha';stageProbs{t-1}.beta]);

if newCutFlag
    stageProbs{t-1}.alpha(end+1,1) = alphan;
    stageProbs{t-1}.beta(:,end+1) = betan;
    stageProbs{t-1}.neg_beta_eta(:,end+1) = sparse([-betan;1]);
else
    fprintf("[L-shaped][It%d] The optimality singlecut existed in Stage-%d: alpha = %.2f, beta = [%s]\n",k,t-1,alphan,sprintf("%.2f ",betan'));
end
ub = stageProbs{t-1}.obj'*x{t-1} + ojbstage2;

out = [];
gap = (ub - lb)/(abs(lb)+1e-6);
if gap < opts.tolerance
    out.stop = 1;
else
    out.stop = 0;
end
fprintf("[L-shaped][It%d] gap = %.2e, lb = %.4f\n",k,gap,lb);
out.ub = ub;
out.lb = lb;
out.x = x{1};
end

function [newCutFlag] = addCut2Pool(coef,mat_exist)
newCutFlag = 1;
m = size(mat_exist,2);
for i = 1:m
    if max(abs(coef(:)-mat_exist(:,i))) < 1e-6
        newCutFlag = 0;
        break
    end
end
end