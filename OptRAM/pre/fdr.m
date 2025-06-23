function newmat=fdr(rawMat,networks)
%networks from 'networks.v7.onlyTF.FDR0.05.xlsx'

[a,b]=size(rawMat.mat);
mat=zeros(a,b);
for i=1:length(networks)
    tf=networks(i,1);
    tar=networks(i,2);
    x=find(strcmp(rawMat.TF,tf),1);
    y=find(strcmp(rawMat.gene,tar),1);
    if ((~isempty(x))&&(~isempty(y)))
        mat(y,x)=1;
    end
end
newmat.mat=rawMat.mat.*mat;
newmat.mat(all(newmat.mat==0,2),:)=[];
newmat.mat(:,all(newmat.mat==0,1))=[];
sumtf=sum(mat,1);
sumgene=sum(mat,2);
j=1;
for i=1:b
    if sumtf(i)~=0
        newmat.TF{j}=rawMat.TF{i};
        j=j+1;
    end
end
j=1;
for i=1:a
    if sumgene(i)~=0
        newmat.gene{j}=rawMat.gene{i};
        j=j+1;
    end
end

end