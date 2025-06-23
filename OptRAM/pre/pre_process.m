function [simplemodel] =  pre_process(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The initial of optimization algorithms
% INPUTS
%model
%
%OUTPUTS
%简化模型
%

[lb,ub,~,~] = fluxVariability(model,0);
t=lb+ub+lb.*ub;%删除的反应
rmMAT=model.S(:,t~=0);
simplemodel.rxns=model.rxns(t~=0);
simplemodel.rev=model.rev(t~=0);
simplemodel.lb=model.lb(t~=0);
simplemodel.ub=model.ub(t~=0);
simplemodel.c=model.c(t~=0);
simplemodel.grRules=model.grRules(t~=0);
simplemodel.rxnNames=model.rxnNames(t~=0);
simplemodel.rxnECNumbers=model.rxnECNumbers(t~=0);
%删除对应反应
%反应对应的无用的代谢物
%反应对应的不调控其他反应的基因
%删除基因对应的调控关系
x=sum(abs(rmMAT),2);%0为删除的代谢物
simplemodel.S=rmMAT(x~=0,:);
simplemodel.mets=model.mets(x~=0);
simplemodel.metCharge=model.metCharge(x~=0);
simplemodel.metNames=model.metNames(x~=0);
simplemodel.metFormulas=model.metFormulas(x~=0);
simplemodel.b=model.b(x~=0);

rgMAT=model.rxnGeneMat(t~=0,:);
y=sum(abs(rgMAT));%0为删除的基因
simplemodel.rxnGeneMat=rgMAT(:,y~=0);
simplemodel.genes=model.genes(y~=0);

end