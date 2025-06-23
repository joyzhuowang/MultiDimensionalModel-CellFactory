function [network,rxnp,genep]=extract_onenet(model,v,regnet,indmax,BPCYids,essentialgene)
%单个改造方案的网络提取
tar=3;
[~,rxn_p,mgene_p,~,~,~,~] = regopt_fitness(indmax,model,v,regnet,BPCYids,tar);
[xgene,~,~,~,~,~] =  check_gene_v4(indmax,model,v,BPCYids,regnet,tar,essentialgene);
[xrxn,~,~,~,~] =  check_rxn_v4(indmax,model,v,regnet,BPCYids,tar);
%network.reg
network.reg={};
allgeneid=find(mgene_p~=1);
genep=[allgeneid,mgene_p(allgeneid);];
n1=numel(allgeneid);
matids=zeros(n1,1);
del=[];
for j=1:n1
    if ~isempty(find(regnet.mgeneid==allgeneid(j), 1))
        matids(j)=find(regnet.mgeneid==allgeneid(j));
    else
        del=[del j];
    end
end
matids(del)=[];

ngene=numel(model.genes);
tf=indmax.gene(indmax.gene>ngene);

if numel(tf)>0
    tf=tf-ngene;
    regnetmat=regnet.mat(matids,tf);
    [row,col]=find(regnetmat~=0);
    n=numel(row);
    fid1=fopen('regnet.txt','w');
    fprintf(fid1,'node1\tnode2\ttype\tcoeff\n');
    for i=1:n
        tfid=tf(col(i));
        mgeneid=matids(row(i));
        regulator=regnet.TFnames{tfid};
        target=regnet.mgene{mgeneid};
        coeff=regnetmat(row(i),col(i));
        network.reg{i,1}=regulator;
        network.reg{i,2}=target;
        if ~isempty(find(xgene==regnet.mgeneid(mgeneid), 1))
            network.reg{i,3}=1;
        else
            network.reg{i,3}=0;
        end
        fprintf(fid1,'%s\t%s\t%d\t%f\n',regulator,target,network.reg{i,3},coeff);
    end
end

%network.met
network.met={};
allrxnid=find(rxn_p~=1);
rxnp=[allrxnid,rxn_p(allrxnid);];
metnetmat=model.rxnGeneMat(allrxnid,allgeneid);
[row,col]=find(metnetmat~=0);
n=numel(row);
fid2=fopen('metnet.txt','w');
fprintf(fid2,'node1\tnode2\ttype\n');
for i=1:n
    rxnid=allrxnid(row(i));
    mgeneid=allgeneid(col(i));
    gene=model.genes{mgeneid};
    rxn=model.rxns{rxnid};
    network.met{i,1}=gene;
    network.met{i,2}=rxn;
    disp(rxnid);
    if ~isempty(find(xrxn==rxnid, 1))
        network.met{i,3}=1;
    else
        network.met{i,3}=0;
    end
    fprintf(fid2,'%s\t%s\t%d\n',gene,rxn,network.met{i,3});
end

fclose('all');
end