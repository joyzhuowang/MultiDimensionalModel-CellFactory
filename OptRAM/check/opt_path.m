function [outs,solution,flux] =  opt_path(indmax,model,v,regnet,BPCYids,index2del,flux,xrxn)

if nargin<=5
    [xrxn,~,~,~,flux] =  check_rxn_v4(indmax,model,v,regnet,BPCYids,3);
end
solution=xrxn;
start=BPCYids(3);
goal=BPCYids(2);
p=flux;
p(abs(p)<0.1)=0;
%outs=1;
q=ones(numel(model.mets),1);
q(index2del)=0;
outs = main( model , start , goal , p , q ,  solution );
end
