function [v,regnet,model] = OptRAM_init(rawmodel,rawregnet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The initialization of optimization algorithms
% INPUTS
%rawmodel--original metabolic model
%rawregnet--original regulatory network
%
% OUTPUTS
% v--reference flux distribution from pFBA 
% model--processed metabolic model
% regnet--processed regulatory network
%

model =  pre_process(rawmodel);
model = pre_calc(model);
%pFBA,1--only minimize the sum of the flux through gene-associated fluxes
[~,~,modelIrrev,MiniFlux] = pFBA(model,'geneoption',1,'skipclass',1);%相当于把rxn中的可逆反应和不可逆反应给区分开来，_f和_b代表可逆反应，_r代表不可逆反应且反向，无标识代表不可逆反应且正向
v=zeros(numel(model.rxns),1);
cnt=1;
for i = 1:length(modelIrrev.rxns)-1
  if (modelIrrev.match(i) ~= 0) %match代表反应_f形式对应的_b形式的位置，例如18号反应的match为1321，说明{'rxn00097_c0_f'}对应的{'rxn00097_c0_b'}为modelirrev.rxns中的1321号反应
    % Reversible reaction
    if (strcmpi(modelIrrev.rxns{i}(end-1:end),'_b'))
      v(cnt) = v(cnt)-MiniFlux.x(i);%可逆反应的v为（-通过该反应的最小flux）
      cnt = cnt + 1;
    else
      v(cnt) = MiniFlux.x(i);
%此处没有给cnt＋1，可以理解为将可逆反应和不可逆反应的v分开，因为rxn00097_c0_f和rxn00097_c0_b等价，所以最后就是不可逆反应的v就是通过该反应的最小flux，可逆反应的v为（-通过该反应的最小flux）
%需要注意的是反应的_f和_b形式的miniflux不一样！
    end
  else
    % Non-reversible reaction，不可逆反应的v就是通过该反应的最小flux
    v(cnt) = MiniFlux.x(i);
    cnt = cnt + 1;
  end
end
%regnet
mat=rawregnet.mat;
targets=rawregnet.mgene;
regnet.TFnames=rawregnet.TFnames;
[newtargets,newmat] = find_mat_mgene(model,mat,targets);
regnet.mat=newmat;
regnet.mgene=newtargets;
tfremove=find(sum(abs(regnet.mat))==0);
regnet.TFnames(tfremove)=[];
regnet.mat(:,tfremove)=[];
generemove=find(sum(abs(regnet.mat),2)==0);
regnet.mgene(generemove)=[];
regnet.mat(generemove,:)=[];
n=numel(regnet.TFnames);
regnet.TFexp=ones(n,1);
regnet.mgeneid=zeros(numel(regnet.mgene),1);
for i=1:numel(regnet.mgene)
    x=find(strcmpi(model.genes,regnet.mgene(i)));
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
    idx=find(strcmpi(targets,model.genes(i)), 1);
    if ~isempty(idx)
        cnt=cnt+1;
        newmat(cnt,:)=mat(idx,:);
    else
        del=[del,i];
    end
end
newtargets(del)=[];
end