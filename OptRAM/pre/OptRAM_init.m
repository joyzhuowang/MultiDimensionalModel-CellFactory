function [v,regnet,model] = OptRAM_init(rawmodel,rawMat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The initial of optimization algorithms
% INPUTS
%rawmodel--原始代谢模型
%rawMat--原始调控网络矩阵
%rawExp--原始表达谱数据
%
% OUTPUTS
% v--初始pFBA得到的参考流值
% model--代谢模型
% regnet--调控网络:matrix;mgene;TFnames;mgeneExp;TFexp
%
model =  pre_process(rawmodel);
model = pre_calc(model);
%pFBA
[~,~,modelIrrev,MiniFlux] = pFBA(model,'geneoption',1,'skipclass',1);
%1--only minimize the sum of the flux through gene-associated fluxes
v=zeros(numel(model.rxns),1);
cnt=1;
for i = 1:length(modelIrrev.rxns)-1
  if (modelIrrev.match(i) ~= 0)
    % Reversible reaction
    if (strcmp(modelIrrev.rxns{i}(end-1:end),'_b'))
      v(cnt) = v(cnt)-MiniFlux.x(i);
      cnt = cnt + 1;
    else
      v(cnt) = MiniFlux.x(i);
    end
  else
    % Non-reversible reaction
    v(cnt) = MiniFlux.x(i);
    cnt = cnt + 1;
  end
end
%regnet读入
mat=rawMat.mat;
targets=rawMat.gene;
regnet.TFnames=rawMat.TF;
[newtargets,newmat] = find_mat_mgene(model,mat,targets);
regnet.mat=newmat;
regnet.mgene=newtargets;
n=numel(regnet.TFnames);
regnet.TFexp=ones(n,1);
%没有加fdr筛选
regnet.mgeneid=zeros(numel(regnet.mgene),1);
for i=1:numel(regnet.mgene)
    x=find(strcmp(model.genes,regnet.mgene(i)));
    regnet.mgeneid(i)=x;
end
end

function [newtargets,newmat] = find_mat_mgene(model,mat,targets)
n=length(model.genes);
newmat=[];
newtargets=model.genes;
cnt=0;
del=[];
for i=1:n
    idx=find(strcmp(targets,model.genes(i)), 1);
    if ~isempty(idx)
        cnt=cnt+1;
        newmat(cnt,:)=mat(idx,:);
    else
        del=[del,i];
    end
end
newtargets(del)=[];
end