%%����һ��ģ���˻�������
%1.�����������, �ļ����ṩ��data�ļ�����
mat=xlsread('CoefMatrix.csv');%���־���,
TF=textread('TF.txt','%s');%�ַ�������
gene=textread('gene.txt','%s');%�ַ�������
[~,networks]=xlsread('fdr005.csv');%�����ַ�������
rawMat.mat=mat;
rawMat.TF=TF;
rawMat.gene=gene;
newmat=fdr(rawMat,networks);
%2.�����л����
initCobraToolbox;
rawmodel=readCbModel('yeast_7.6_cobra.xml');%��лģ��
%3.���������ʼ��
[v,regnet,model] = OptRAM_init(rawmodel,newmat);
%4.BPCYids����
BPCYids=zeros(3,1);
BPCYids(1)=find(strcmp(model.rxnNames,'growth'));%ϸ��������Ӧ��,�������Ŀ������
BPCYids(2)=find(strcmp(model.rxnNames,'succinate exchange'));%Ŀ����ﷴӦ��,�������Ŀ������
BPCYids(3)=find(strcmp(model.rxnNames,'D-glucose exchange'));%Ӫ�����ﷴӦ��,�������Ŀ������
%5.������
save data.mat model v regnet BPCYids
[new_model,fitmax,indmax] =  OptRAM_main(model,v,regnet,BPCYids,'YOURPATH');%���logָ��λ��