function [fitness,rxn_p,mgene_p,growth,target,model2] =  Optfitness_mgene(individual,model,v,BPCYids)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%rxn_p,model,FBAsolution,growth,target,substrate
% The fitness function of optimized individuals
% INPUTS individual,model,v,regulations 
%
%
mgene_p=ones(numel(model.genes),1);
rules=cell(numel(model.grRules),1);
for i=1:numel(model.grRules)
    rules{i}=model.grRules{i};
end
expression=model.postfix_expression;
long=model.postfix_long;
%从individual的code到mgene
for i=1:length(individual.gene)
    mgene_p(individual.gene(i))=precode(individual.code(i));
end

%从gene的code到rxn
n=length(model.rxns);
rxn_p=ones(n,1);
%对全部反应计算其p值，优化目标：只计算改变的反应
%for i=1:n
%    if ~isempty(model.grRules{i})
%        rxn_p(i) = calc( model,model.grRules{i},mgene_p);
%    end
%end
for i=1:n
    if ~isempty(rules{i})
        temp = calc_postfix(expression,long,i,mgene_p);
        rxn_p(i) = temp(1);
    end
end
%modelWT=model;
%根据rxn的p改变约束
%rxn_p从上面继承，为二维数组，由改变的rxnid和p组成
%model.lb=lb;
%model.ub=ub;
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
model = changeObjective(model,model.rxns(BPCYids(1)));
FBAsolution = optimizeCbModel(model);%默认目标反应为biomass
%[FBAsolution,~,~,~] = MOMA(modelWT,model);
if numel(FBAsolution.x)==0
    fitness=0;
    growth=0;
    target=0;
else
    model.lb(BPCYids(1)) = FBAsolution.f*0.99;
    model2 = changeObjective(model,model.rxns(BPCYids(2)));
            FBAsolution1 = optimizeCbModel(model2,'min');
            FBAsolution2 = optimizeCbModel(model2,'max');
            if numel(FBAsolution1.x)==numel(FBAsolution2.x) && numel(FBAsolution1.x)~=0
                FBAsolution.x = (FBAsolution1.x+FBAsolution2.x)./2;
                rang = (FBAsolution2.x(BPCYids(2))-FBAsolution1.x(BPCYids(2)))/2;
            else
                FBAsolution.x = [];
                rang=100;
            end   
    if numel(FBAsolution.x)==0
        fitness=0;
        growth=0;
        target=0;
    else
        growth=FBAsolution.x(BPCYids(1));%growthid
        target=FBAsolution.x(BPCYids(2));%targetid
        substrate=FBAsolution.x(BPCYids(3));%substrateid
        fitness=abs(growth*target/substrate)*(1-log(rang/target));
        %x1 fitness=abs(growth*target/substrate)/(range/target);
        %fitness=abs(growth*target/substrate);
    end
end

end

function p=precode(code)
if code==0
    p=0.001;
elseif code==6
    p=1;
else
    p=2^code;
end
end
