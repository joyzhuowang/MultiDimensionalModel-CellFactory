%���������������
%load('data.mat')

%index2del
load('smallmol.mat')
%essentialgene
load('essential.mat')
%����log�ļ�����
nlog=10;

[summary,effects,xgenes,xrxns,geneids,codes,genes,counttable,fluxes,paths,pathoutput]=multi_check(model,v,regnet,BPCYids,index2del,essentialgene,nlog);

file=('output.xlsx');
%���ܱ�
xlswrite(file,summary,1)
%����λ��
xlswrite(file,genes,2);
xlswrite(file,codes,3);
%��ֱ�
xlswrite(file,effects,4)
%·�����
save pathresult.mat pathoutput

%����ؼ�λ�㹹�ɵ�����
[counttable]=extract_allnet(model,regnet,counttable,xrxns,xgenes);