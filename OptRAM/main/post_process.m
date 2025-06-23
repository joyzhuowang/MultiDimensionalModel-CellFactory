function [simpleindmax,fit] =  post_process(indmax,new_model,v,regnet,BPCYids,tar)
%����ȥ����Ӱ���ֵ�λ��
% INPUTS
%indmax--ԭʼ�ĸ������
%model--��лģ��
%v--�ο���ֵ
%regnet--��������:matrix;mgene;TFnames;mgeneExp;TFexp
%BPCYids--BPCY�з�Ӧ��ID
%tar--Ŀ�꺯��ѡ�����
%
% OUTPUTS
% simpleindmax--�򻯺�ĸ������
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