function [xgene,xgene_p,result,final_fit,flux,essdown] =  check_genes(indmax,model,v,BPCYids,regnet,tar,essentialgene)
if nargin>4
    [~,~,mgene_p,~,~,~,~] = OptRAM_obj(indmax,model,v,regnet,BPCYids,tar);
else
    [~,~,mgene_p,~,~] =  Optfitness_mgene(indmax,model,v,BPCYids);    
end
xgene=find(mgene_p~=1);
xgene_p=mgene_p(xgene);
essdown=0;
[matchgene,id]=intersect(xgene,essentialgene);
if ~isempty(matchgene)
    for j=1:numel(matchgene)
        if xgene_p(id(j))<0.5
            essdown=essdown+1;
        end
    end
end
[fit,target,growth,flux,~,mintar,maxtar] = fitness(model,v,BPCYids,mgene_p,tar);
final_fit=[fit,target,growth,mintar,maxtar];
x1=round(fit*10000)/10000;
n=numel(xgene);
result=zeros(n,6);
result(:,1)=xgene;
result(:,2)=xgene_p;
delid=[];
for i=1:n
    temp_mgene_p=mgene_p;
    temp_mgene_p(xgene(i))=1;
    [fit,~,growth,~,~,mintar,maxtar] = fitness(model,v,BPCYids,temp_mgene_p,tar);
    result(i,3)=fit/x1;
    result(i,4)=mintar;
    result(i,5)=maxtar;
    result(i,6)=growth;
    x2=round(fit*10000)/10000;
    if x2>=x1*0.9
        delid=[delid i];
    end
end 
result(delid,:)=[];
xgene(delid)=[];
xgene_p(delid)=[];
%disp(find(new_mgene_p~=1));
%disp(fit);
end


function [fit,target,growth,flux,range,mintar,maxtar] = fitness(model,v,BPCYids,mgene_p,tar)
rules=model.grRules;
expression=model.postfix_expression;
long=model.postfix_long;
n=length(model.rxns);
rxn_p=ones(n,1);
for i=1:n
    if ~isempty(rules{i})
        temp = calc_postfix(expression,long,i,mgene_p);
        rxn_p(i) = temp(1);
    end
end
for i=1:n
    p=rxn_p(i);
    if p>1
        if v(i)>0
            model.lb(i)=min(p*v(i),model.ub(i));
        elseif v(i)<0
            model.ub(i)=max(p*v(i),model.lb(i));
        end
    elseif p<1
        if v(i)>0
            model.ub(i)=p*v(i);
            model.lb(i)=0;
        elseif v(i)<0
            model.lb(i)=p*v(i);
            model.ub(i)=0;
        else
            model.lb(i)=0;
            model.ub(i)=0;
        end
    end
end

FBAsolution = optimizeCbModel(model);
%[FBAsolution,~,~,~] = MOMA(modelWT,model);
flux=0;
range=0;
mintar=0;
maxtar=0;
if  numel(FBAsolution.x)==0
    fit=0;
    growth=0;
    target=0;
else
    model.lb(BPCYids(1)) = FBAsolution.f*0.99;
    model2 = changeObjective(model,model.rxns(BPCYids(2)));
    range=0;
    switch tar
        case 1
            FBAsolution = optimizeCbModel(model2,'min');
        case 2
            FBAsolution = optimizeCbModel(model2,'max');
        case 3
            FBAsolution1 = optimizeCbModel(model2,'min');            
            FBAsolution2 = optimizeCbModel(model2,'max');         
            if numel(FBAsolution1.x)==numel(FBAsolution2.x) && numel(FBAsolution1.x)~=0
                FBAsolution.x = (FBAsolution1.x+FBAsolution2.x)./2;
                range=(FBAsolution2.x(BPCYids(2))-FBAsolution1.x(BPCYids(2)))/2;
            else
                FBAsolution.x = [];
                range=100;
            end
    end 
    if numel(FBAsolution.x)==0
        fit=0;
        growth=0;
        target=0;
        flux=0;
    else
        growth=FBAsolution.x(BPCYids(1));%growthid
        target=FBAsolution.x(BPCYids(2));%targetid
        substrate=FBAsolution.x(BPCYids(3));%substrateid
        fit=abs(growth*target/substrate)*(1-log(range/target));
        flux=FBAsolution.x;
        mintar=FBAsolution1.x(BPCYids(2));
        maxtar=FBAsolution2.x(BPCYids(2));
    end
end

end