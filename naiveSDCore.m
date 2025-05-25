function [stageProbs,bundlePI,bundleJ,OmegaInfo,out] = ...
        naiveSDCore(stageProbs,rpm,OmegaInfo,bundlePI,bundleJ,observ,k,constL,opts,x1)
% Stochastic decomposition without regularization
% Inputs:
% stageProbs{t}: the main program at stage t
% rpm: uncertainty at stage t
% OmegaInfo: sampling information at stage 2
% bundlePI: the explored dual solution at stage 2 
% bundleJ: the cuts approximating cost-to-go function
% observ: new obervation
% k: iteration index
% constL: the lower bound of cost-to-go function
% opts: algorithm options
% - opts.verbose: 1 for printing the iteration information
% - opts.logdir: iterate log save path
% x1: initial iterate at stage 1
% 
% Outputs:
% out.ub: the upper bound at iteration k
% out.lb: the lower bound at iteration k
% out.x: the solution at iteration k
% 
% Author: Zhiyuan Zhang
% Date: May 25, 2025


stageNum = length(stageProbs);
params = [];
params.outputflag = 0;
params.InfUnbdinfo = 1;
x = cell(stageNum,1);

% Forward Pass along with the scenario path
if opts.verbose
    fprintf("[It%d] observ = %s\n",k,sprintf("%d ",observ));
end

%% Root-stage
t = 1;
[m1,n1] = size(stageProbs{t}.D);

if ~isempty(x1)
    x{t} = x1;
    lb = stageProbs{t}.obj'*x{t} + constL;
else
    numeta = length(OmegaInfo.indices);
    programMain = [];
    if isfield(stageProbs{t},'ub')
        programMain.ub = stageProbs{t}.ub;
        if numeta >= 1
            programMain.ub(end+1,1) = inf;
        end
    end

    if isfield(stageProbs{t},'lb')
        programMain.lb = stageProbs{t}.lb;
        if numeta >= 1
            programMain.lb(end+1,1) = constL;
        end
    end

    if numeta >= 1
        programMain.obj = [stageProbs{t}.obj; 1];
        programMain.varnames = stageProbs{t}.varnames;
        programMain.varnames{end+1} = 'eta';

        % scale the bundle set
        coef = 1 - 1/k;
        bundleJ.alpha(1:end-1) = coef * bundleJ.alpha(1:end-1) + (1-coef) * constL;
        bundleJ.beta(1:end-1) = coef * bundleJ.beta(1:end-1);

        rhs_add = bundleJ.alpha;
        At_add = [-bundleJ.beta',ones(bundleJ.cnt,1)];

        A_r_add = unique([At_add,rhs_add], 'rows');
        rhs_add = A_r_add(:,end);
        At_add = A_r_add(:,1:end-1);
        programMain.A = [stageProbs{t}.D, zeros(m1,1); At_add];
        programMain.rhs = [stageProbs{t}.r; rhs_add];
        sensecnt = size(At_add,1);
    else
        programMain.varnames = stageProbs{t}.varnames;
        programMain.obj = stageProbs{t}.obj;
        programMain.A = stageProbs{t}.D;
        programMain.rhs = stageProbs{t}.r;
        programMain.objcon = constL; % for the cost-to-go function value
        sensecnt = 0;
    end

    programMain.sense = strcat(stageProbs{t}.sense,repmat('>',1,sensecnt));


    fileName = fullfile(opts.logdir,sprintf('It%dFPstage%ddProb.lp',k,t));
    if opts.writeflag
        gurobi_write(programMain, fileName);
    end
    outxt = gurobi(programMain, params);

    if strcmp(outxt.status,'OPTIMAL')
        x{t} = outxt.x(1:n1);
        lb = outxt.objval;
    end
end

if opts.verbose
    fprintf("[It%d] lb = %.6f.\n", k, lb);
end


% add new sample to Omega
newOmegaFlag = 1;
for j = 1:length(OmegaInfo.indices)
    if observ == OmegaInfo.indices(j)
        newOmegaFlag = 0;
        OmegaInfo.weights(j,1) = OmegaInfo.weights(j,1) + 1;
        break
    end
end

if newOmegaFlag
    OmegaInfo.indices(end+1) = observ;
    OmegaInfo.weights(end+1,1) = 1;
end
OmegaInfo.cnt = OmegaInfo.cnt + 1;


% Backward recursion
t = 2;

programMain = [];
programMain.A = stageProbs{t}.D;
programMain.sense = stageProbs{t}.sense;
programMain.varnames = stageProbs{t}.varnames;
programMain.obj = stageProbs{t}.obj;
rti_observ = stageProbs{t}.r + rpm.val(:,observ);
programMain.rhs = rti_observ - stageProbs{t}.C * x{t-1};

if isfield(stageProbs{t},'ub') && isfield(stageProbs{t},'lb')
    programMain.lb = stageProbs{t}.lb;
    programMain.ub = stageProbs{t}.ub;
elseif isfield(stageProbs{t},'ub')
    programMain.ub = stageProbs{t}.ub;
elseif isfield(stageProbs{t},'lb')
    programMain.lb = stageProbs{t}.lb;
end


outxti = gurobi(programMain, params);

if strcmp(outxti.status,'UNBOUNDED')
    warning("subproblem is unbounded!")
end


if strcmp(outxti.status,'OPTIMAL')
    dualvar = [];
    dualvar.pi = outxti.pi;
    if isfield(programMain,'ub') && isfield(programMain,'lb')
        idxub = programMain.ub < inf;
        idxlb = programMain.lb > -inf;
        dualvar.dualbd = programMain.ub(idxub)'*min(outxti.rc(idxub),0) + programMain.lb(idxlb)'*max(outxti.rc(idxlb),0);
    elseif isfield(programMain,'ub')
        idxub = programMain.ub < inf;
        dualvar.dualbd = programMain.ub(idxub)'*min(outxti.rc(idxub),0);
    elseif isfield(programMain,'lb')
        idxlb = programMain.lb > -inf;
        dualvar.dualbd = programMain.lb(idxlb)'*max(outxti.rc(idxlb),0);
    else
        dualvar.dualbd = 0;
    end

    num_sce = length(OmegaInfo.indices);
    tmp_alpha = zeros(num_sce,1);
    tmp_beta = zeros(n1,num_sce);
    tmp_obj = zeros(num_sce,1);
    bundlePI = updateBundleSet(k,dualvar,bundlePI);
    Cx_1 = stageProbs{t}.C * x{t-1};
    for i = 1: num_sce
        idx = OmegaInfo.indices(i);
        if idx == observ
            tmp_alpha(i,1) = rti_observ' * outxti.pi + dualvar.dualbd;
            tmp_beta(:,i) = - stageProbs{t}.C' * outxti.pi;
            tmp_obj(i,1) = outxti.objval;
        else
            % argmax
            rti = stageProbs{t}.r + rpm.val(:,idx);
            rhs = rti - Cx_1;
            [valmax,idxmax] = max(bundlePI.pi' * rhs + bundlePI.dualbd);
            pimax = bundlePI.pi(:,idxmax);
            tmp_alpha(i,1) = rti' * pimax + bundlePI.dualbd(idxmax,1);
            tmp_beta(:,i) = - stageProbs{t}.C' * pimax;
            tmp_obj(i,1) = valmax;
        end
    end

    probvec = OmegaInfo.weights / OmegaInfo.cnt;

    ub = tmp_obj' * probvec;
    alphan = tmp_alpha' * probvec;
    betan = tmp_beta * probvec;

    newCutFlag = 1;
    if newCutFlag
        bundleJ.alpha(end+1,1) = alphan;
        bundleJ.beta(:,end+1) = betan;
        bundleJ.cnt = bundleJ.cnt + 1;
    end

end

out = [];
if opts.verbose
    fprintf("[It%d] Complete the backward recursion.\n%s\n",k,repmat('*',1,20));
    fprintf("[It%d] ub = %.6f.\n", k, ub);
    gap = ub - lb;
    fprintf("[It%d] ub - lb = %.2e, (ub+lb)/2 = %.6e\n",k,gap,0.5*(ub+lb));
end
out.ub = ub;
out.lb = lb;
out.x = x{1};
end


% function [status] = PrintSolution(fileName,model,out)
% fileID = fopen(fileName, 'a');
% fprintf(fileID, 'Solution\n');
% for j = 1 : size(model.A,2)
%     fprintf(fileID, '%s = %f\n',model.varnames{j},out.x(j));
% end
% status = fclose(fileID);
% end

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


function [bundlePI] = updateBundleSet(k,vpi,bundlePI)
flag = 1;
[~,m] = size(bundlePI.pi);
for j = 1:m
    
    if max(abs([vpi.pi - bundlePI.pi(:,j);vpi.dualbd - bundlePI.dualbd(j,1)])) < 1e-6
        flag = 0;
        break
    end
end
if flag == 1
    bundlePI.pi(:,end+1) = vpi.pi;
    bundlePI.dualbd(end+1,1) = vpi.dualbd;
    bundlePI.size = bundlePI.size + 1;
    bundlePI.iter(end+1,1) = k;
end
end

function [output,flag] = updateIndexSet(idx_set,new_idx)
output = idx_set;
[~,m] = size(idx_set);
flag = 1;
for j = 1:m
    if norm(new_idx - idx_set(:,j)) < 1e-6
        flag = 0;
        break
    end
end
if flag == 1
    output(:,end+1) = new_idx;
end
end