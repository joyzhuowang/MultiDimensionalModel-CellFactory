%后续处理，结果汇总
%load('data.mat')

%index2del
load('smallmol.mat')
%essentialgene
load('essential.mat')
%输入log文件个数
nlog=10;

[summary,effects,xgenes,xrxns,geneids,codes,genes,counttable,fluxes,paths,pathoutput]=multi_check(model,v,regnet,BPCYids,index2del,essentialgene,nlog);

file=('output.xlsx');
%汇总表
xlswrite(file,summary,1)
%具体位点
xlswrite(file,genes,2);
xlswrite(file,codes,3);
%打分表
xlswrite(file,effects,4)
%路径输出
save pathresult.mat pathoutput

%输出关键位点构成的网络
[counttable]=extract_allnet(model,regnet,counttable,xrxns,xgenes);