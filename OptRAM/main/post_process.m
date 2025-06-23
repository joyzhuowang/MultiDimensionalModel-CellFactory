function [simpleindmax,fit] =  post_process(indmax,new_model,v,regnet,BPCYids,tar)
%后处理，去除不影响打分的位点
% INPUTS
%indmax--原始的改造个体
%model--代谢模型
%v--参考流值
%regnet--调控网络:matrix;mgene;TFnames;mgeneExp;TFexp
%BPCYids--BPCY中反应的ID
%tar--目标函数选择参数
%
% OUTPUTS
% simpleindmax--简化后的改造个体
%

if numel(indmax.gene)==0
    simpleindmax=0;
    fit=0;
else
[fitmax,~,~,target,growth,~,~]=regopt_fitness(indmax,new_model,v,regnet,BPCYids,tar);
x1=round(fitmax*10000)/10000;
n=numel(indmax.gene);
fit=zeros(n+1,3);
delid=[];
fit(1,1)=x1;
fit(1,2)=target;
fit(1,3)=growth;
for i=1:n
    tempindmax=indmax;
    tempindmax.code(i)=6;
    [tempfit,~,~,target,growth]=regopt_fitness(tempindmax,new_model,v,regnet,BPCYids,tar);
    fit(i+1,1)=tempfit;
    fit(i+1,2)=target;
    fit(i+1,3)=growth;
    x2=round(tempfit*10000)/10000;
    if x2>=x1*0.99
        delid=[delid i];
    end
end 
simpleindmax=indmax;
simpleindmax.gene(delid)=[];
simpleindmax.code(delid)=[];
fit(delid,:)=[];
end
end