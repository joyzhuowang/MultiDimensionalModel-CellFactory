function [simplemodel] =  pre_process(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The initial of optimization algorithms
% INPUTS
%model
%
%OUTPUTS
%��ģ��
%

[lb,ub,~,~] = fluxVariability(model,0);
t=lb+ub+lb.*ub;%ɾ���ķ�Ӧ
rmMAT=model.S(:,t~=0);
simplemodel.rxns=model.rxns(t~=0);
simplemodel.rev=model.rev(t~=0);
simplemodel.lb=model.lb(t~=0);
simplemodel.ub=model.ub(t~=0);
simplemodel.c=model.c(t~=0);
simplemodel.grRules=model.grRules(t~=0);
simplemodel.rxnNames=model.rxnNames(t~=0);
simplemodel.rxnECNumbers=model.rxnECNumbers(t~=0);
%ɾ����Ӧ��Ӧ
%��Ӧ��Ӧ�����õĴ�л��
%��Ӧ��Ӧ�Ĳ�����������Ӧ�Ļ���
%ɾ�������Ӧ�ĵ��ع�ϵ
x=sum(abs(rmMAT),2);%0Ϊɾ���Ĵ�л��
simplemodel.S=rmMAT(x~=0,:);
simplemodel.mets=model.mets(x~=0);
simplemodel.metCharge=model.metCharge(x~=0);
simplemodel.metNames=model.metNames(x~=0);
simplemodel.metFormulas=model.metFormulas(x~=0);
simplemodel.b=model.b(x~=0);

rgMAT=model.rxnGeneMat(t~=0,:);
y=sum(abs(rgMAT));%0Ϊɾ���Ļ���
simplemodel.rxnGeneMat=rgMAT(:,y~=0);
simplemodel.genes=model.genes(y~=0);

end