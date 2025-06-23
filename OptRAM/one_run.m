%%运行一次模拟退火主函数
%1.读入调控网络, 文件已提供在data文件夹中
mat=xlsread('CoefMatrix.csv');%数字矩阵,
TF=textread('TF.txt','%s');%字符串数组
gene=textread('gene.txt','%s');%字符串数组
[~,networks]=xlsread('fdr005.csv');%两列字符串数组
rawMat.mat=mat;
rawMat.TF=TF;
rawMat.gene=gene;
newmat=fdr(rawMat,networks);
%2.读入代谢网络
initCobraToolbox;
rawmodel=readCbModel('yeast_7.6_cobra.xml');%代谢模型
%3.整合网络初始化
[v,regnet,model] = OptRAM_init(rawmodel,newmat);
%4.BPCYids设置
BPCYids=zeros(3,1);
BPCYids(1)=find(strcmp(model.rxnNames,'growth'));%细胞生长反应名,根据你的目标来改
BPCYids(2)=find(strcmp(model.rxnNames,'succinate exchange'));%目标产物反应名,根据你的目标来改
BPCYids(3)=find(strcmp(model.rxnNames,'D-glucose exchange'));%营养底物反应名,根据你的目标来改
%5.主函数
save data.mat model v regnet BPCYids
[new_model,fitmax,indmax] =  OptRAM_main(model,v,regnet,BPCYids,'YOURPATH');%输出log指定位置