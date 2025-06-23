function [indmax,xgene,xrxn,result1,result2,summary,effect,flux] =  check_result(file,model,v,regnet,BPCYids,essentialgene)
%indmax,model,v,regnet,BPCYids
%%自动读入每个个体
%%对个体分析，至有效基因，有效反应，途径寻找
model_pre=pre_calc(model);

fidin=fopen(file,'r');
fseek(fidin,-100,'eof');
fgetl(fidin);
while ~feof(fidin)
    x=fgetl(fidin);
end
y=x;
s=regexp(y,'\t','split');
ind_num=(numel(s)-3)/2;
indmax.gene=zeros(1,ind_num);
indmax.code=zeros(1,ind_num);
for i=1:ind_num
    indmax.gene(i)=str2num(s{i});
end
for i=1:ind_num
    indmax.code(i)=str2num(s{ind_num+i});
end
[xgene,~,result1,fit2,~,essdown] =  check_genes(indmax,model_pre,v,BPCYids,regnet,3,essentialgene);
%[table] = gene2tf(xgene,indmax,regnet);
[xrxn,~,result2,~,flux] =  check_rxns(indmax,model_pre,v,regnet,BPCYids,3);
%x=model.rxnGeneMat(xrxn,xgene);
%mat=full(x);

variation=flux'*v/(norm(flux)*norm(v));
%variation=norm(flux-v);
summary=zeros(1,7);
summary(1,1)=fit2(1);
summary(1,2)=fit2(4);
summary(1,3)=fit2(5);
summary(1,4)=fit2(3);
summary(1,5)=numel(find(indmax.gene<=714));
summary(1,6)=numel(find(indmax.gene>714));
summary(1,7)=variation;

effect=[numel(xgene),numel(xrxn),essdown];

end