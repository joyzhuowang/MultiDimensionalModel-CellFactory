function [summary,effects,xgenes,xrxns,geneids,codes,genes,counttable,fluxes,paths,pathoutput]=multi_check(model,v,regnet,BPCYids,index2del,essentialgene,nlog)
%�����xgenes��xrxns��Ϊ����������Ӱ�쳬��10%�ģ����һ��ΪӰ���ɾȥ��λ�����������Ϊ���ŷ����İٷֱ�
summary=zeros(nlog,7);
effects=zeros(nlog,3);
paths=zeros(nlog,5);
fluxes=[];
xrxns=[];
xrxnsr=[];
xgenes=[];
xgenesr=[];
genes={};
geneids=zeros(nlog,10);
codes=zeros(nlog,10);
pathoutput={};
for i=1:nlog
    filename=['log',num2str(i),'.txt'];
    disp(filename);
    [indmax,~,xrxn,result1,result2,sum,effect,flux] =  check_result(filename,model,v,regnet,BPCYids,essentialgene);
    [outs,solution,~] =  opt_path(indmax,model,v,regnet,BPCYids,index2del,flux,xrxn);
    pathl=outs.path_length;
    pathf=outs.path_flow;
    branchm=mean(outs.branch_length);
    score=outs.pathScore;
    pathoutput{i}=outs;
    n_sol=numel(solution);
    n=numel(indmax.gene);
    for j=1:n
        if indmax.gene(j)<=numel(model.genes)
            genes{i,j}=model.genes(indmax.gene(j));
        else
            genes{i,j}=regnet.TFnames(indmax.gene(j)-numel(model.genes));
        end
    end
    geneids(i,1:n)=indmax.gene;
    codes(i,1:n)=indmax.code;
    xgenes=[xgenes result1(:,1)' 0];
    xgenesr=[xgenesr result1(:,3)' 0];
    xrxns=[xrxns result2(:,1)' 0];
    xrxnsr=[xrxnsr result2(:,3)' 0];
    summary(i,:)=sum;
    effects(i,:)=effect;
    fluxes(:,i)=flux;
    paths(i,1)=pathl;
    paths(i,2)=pathf;
    paths(i,3)=branchm;
    paths(i,4)=n_sol;
    paths(i,5)=score;
end
xrxns=[xrxns;xrxnsr];
xrxns=xrxns';
xgenes=[xgenes;xgenesr];
xgenes=xgenes';
genevect=geneids(:);
n_allgene=numel(model.genes)+numel(regnet.TFnames);
genecount=histc(genevect,1:n_allgene);
n=numel(find(genecount>0));
counttable=zeros(n,2);
counttable(:,1)=find(genecount>0);
counttable(:,2)=genecount(genecount>0);
end