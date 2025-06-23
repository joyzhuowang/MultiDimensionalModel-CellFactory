function [xrxn,xrxn_p,result,final_fit,flux] =  check_rxns(indmax,model,v,regnet,BPCYids,tar)
if nargin>4
    [~,rxn_p,~,~,~,~,model2] = OptRAM_obj(indmax,model,v,regnet,BPCYids,tar);
else
    [~,rxn_p,~,~,~] =  Optfitness_mgene(indmax,model,v,BPCYids);    
end
%[~,~,~,fluxmax] = fitness(model,v,BPCYids,rxn_p);
%get flux with pFBA
model2.c(BPCYids(1))=1;
model2.lb(BPCYids(1))=0;
model2.c(BPCYids(2))=0;
[~,~,modelIrrev,MiniFlux] = pFBA(model2,'geneoption',1,'skipclass',1);
x=zeros(numel(model2.rxns),1);
cnt=1;
if numel(MiniFlux.x)~=0
    for j = 1:length(modelIrrev.rxns)-1
        if (modelIrrev.match(j) ~= 0)
            % Reversible reaction
            if (strcmp(modelIrrev.rxns{j}(end-1:end),'_b'))
                x(cnt) = x(cnt)-MiniFlux.x(j);
                cnt = cnt + 1;
            else
                x(cnt) = MiniFlux.x(j);
            end
        else
            % Non-reversible reaction
            x(cnt) = MiniFlux.x(j);
            cnt = cnt + 1;
        end
    end
else
    disp('warn pFBA');
end
flux=x;

xrxn=find(rxn_p~=1);
xrxn_p=rxn_p(xrxn);
[fit,~,growth,~,~,mintar,maxtar] = fitness(model,v,BPCYids,rxn_p,tar);
x1=round(fit*10000)/10000;
final_fit=[x1,growth,mintar,maxtar];
n=numel(xrxn);
result=zeros(n,6);
result(:,1)=xrxn;
result(:,2)=xrxn_p;
delid=[];
delid2=[];
for i=1:n
    temp_rxn_p=rxn_p;
    temp_rxn_p(xrxn(i))=1;
    [fit,~,growth,~,~,mintar,maxtar] = fitness(model,v,BPCYids,temp_rxn_p,tar);
    result(i,3)=fit/x1;
    result(i,4)=mintar;
    result(i,5)=maxtar;
    result(i,6)=growth;
    x2=round(fit*10000)/10000;
    if x2>=x1*0.9;
        delid=[delid i];
    end
    if x2>=x1*0.9;
        delid2=[delid2 i];
    end
end
result(delid,:)=[];
xrxn(delid2)=[];
xrxn_p(delid2)=[];

%disp(x1);
%disp(numel(xrxn));
end


function [fit,target,growth,flux,range,mintar,maxtar] = fitness(model,v1,BPCYids,rxn_p,tar)
n=numel(rxn_p);
for i=1:n
    p=rxn_p(i);
    if p>1
        if v1(i)>0
            model.lb(i)=min(p*v1(i),model.ub(i));
        elseif v1(i)<0
            model.ub(i)=max(p*v1(i),model.lb(i));
        end
    elseif p<1
        if v1(i)>0
            model.ub(i)=p*v1(i);
            model.lb(i)=0;
        elseif v1(i)<0
            model.lb(i)=p*v1(i);
            model.ub(i)=0;
        else
            model.lb(i)=0;
            model.ub(i)=0;
        end
    end
end
flux=0;
range=0;
mintar=0;
maxtar=0;
FBAsolution = optimizeCbModel(model);
%[FBAsolution,~,~,~] = MOMA(modelWT,model);
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
