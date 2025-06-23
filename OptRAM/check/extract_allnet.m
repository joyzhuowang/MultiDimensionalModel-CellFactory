function [counttable]=extract_allnet(model,regnet,counttable,xrxns,xgenes)
%整体方案的网络提取
xrxns=xrxns(:,1);
xgenes=xgenes(:,1);
tmpgene=unique(xgenes);
gids=tmpgene(tmpgene>0);
n1=numel(gids);
ids=zeros(n1,1);
del=[];
for j=1:n1
    if ~isempty(find(regnet.mgeneid==gids(j), 1))
        ids(j)=find(regnet.mgeneid==gids(j));
    else
        del=[del j];
    end
end
ids(del)=[];
gids(del)=[];
%找出所有的TF
tf=counttable(counttable(:,1)>714,1);
tf=tf-714;
netmat=regnet.mat(ids,tf);
[row,col]=find(netmat~=0);
n=numel(row);
fid=fopen('net.txt','w');
fprintf(fid,'node1\tnode2\ttype\n');
for i=1:n
    tfid=tf(col(i));
    mgeneid=ids(row(i));
    regulator=regnet.TFnames{tfid};
    target=regnet.mgene{mgeneid};
    coeff=netmat(row(i),col(i));
    fprintf(fid,'%s\t%s\t%f\n',regulator,target,coeff);
end
%mgenes=find(sum(netmat,2)>0);
%mids=regnet.mgeneid(mgenes);
mids=tmpgene(tmpgene>0);
tmprxn=unique(xrxns);
rids=tmprxn(tmprxn>0);
%disp(rids);
netmat2=model.rxnGeneMat(rids,mids);
[row2,col2]=find(netmat2~=0);
n2=numel(row2);
for i=1:n2
    mgeneid=mids(col2(i));
    rxnid=rids(row2(i));
    gene=model.genes{mgeneid};
    rxn=model.rxns{rxnid};    
    fprintf(fid,'%s\t%s\t0\n',gene,rxn);
end

fclose('all');
end