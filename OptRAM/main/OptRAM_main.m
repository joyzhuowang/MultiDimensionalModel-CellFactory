function [new_model,fitmax,indmax] =  OptRAM_main(model,v,regnet,BPCYids,loc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The main optimization algorithm with simulated annealing
% INPUTS
%model--代谢模型
%v--参考流值
%regnet--调控网络:matrix;mgene;TFnames;mgeneExp;TFexp
%BPCYids--BPCY中反应的ID
%tar--目标函数选择参数
%
% OUTPUTS
% fitmax--最大目标函数
% indmax--最大目标函数对应改造个体
%
tar=3;%tar表示target设为目标最大或最小
%初始化模拟退火参数
k=5;
Eo=k/1000; 
%Ef=k/100000;
T_max=-Eo/log(0.5); 
%T_min=-Ef/log(0.5);
iter_max=500; 
%NFEs=400000;%13hours
%alpha=exp((log(T_min)-log(T_max))/(NFEs/iter_max));
alpha=0.995;
s_max=100;
T=T_max;
initial_size=4;
ctime=datestr(now,30);
tseed=str2num(ctime((end-5):end));%设置随机种子
rand('seed',tseed)%主机专用
ntf=numel(regnet.TFnames);
n_mgene=numel(model.genes);
ngene=ntf+n_mgene;
individual1=rand_ind(ngene,ntf,initial_size);%随机生成一个个体的函数
individual11=checkTF(individual1,regnet,n_mgene);
fit1=OptRAM_obj(individual11,model,v,regnet,BPCYids,tar);
filename=[loc,'log.txt'];
fid=fopen(filename, 'w');
fitmax=0;
indmax.gene=0;
indmax.code=0;
flag_num=0;%若解的分数始终未上涨，标记该温度,并记录连续多少个温度未上涨
flag=0;
disp(ntf);
%while flag_num<100 || (fitmax==0 && flag_num<300)
while flag_num<15 || (fitmax==0 && flag_num<30)
    tic
    iter_num=1; 
    s_num=1;
    fprintf(fid,'T:%.5e\n',T);
    while iter_num<iter_max && s_num<s_max; 
        individual2=mutant(individual1,ngene,regnet);%个体突变函数
        individual22=checkTF(individual2,regnet,n_mgene);
        size=numel(individual22.gene);
        %disp(size);
        individual3=rand_ind(ngene,ntf,size);
        individual33=checkTF(individual3,regnet,n_mgene);
        fit2=OptRAM_obj(individual22,model,v,regnet,BPCYids,tar);
        fit3=OptRAM_obj(individual33,model,v,regnet,BPCYids,tar);
        if fit3>=fit2
            individual22=individual33;
            fit2=fit3;
        end
        R=rand;         
        Deltafit=fit2-fit1;
        individual2=individual22;
        %disp(Deltafit);
        if exp(Deltafit/T)>R;%突变后个体分数大，则接受突变
            individual1=individual2; 
            fit1=fit2;
            x1=round(fit1*1000000)/1000000;%判断最大值的计算精度设定为小数点后6位
            s_num=1;
            if x1>fitmax%记录分数最大size最小的个体
                fitmax=x1;
                indmax=individual1;
                for z=1:numel(individual1.gene)
                    fprintf(fid,'%d(%d)\t',individual1.gene(z),individual1.code(z));
                end
                fprintf(fid,'%f\n',fit1);
                flag=1;
            elseif x1==fitmax
                if numel(indmax.gene)>numel(individual1.gene)
                    indmax=individual1;
                    for z=1:numel(individual1.gene)
                        fprintf(fid,'%d(%d)\t',individual1.gene(z),individual1.code(z)); 
                    end
                    fprintf(fid,'%f\n',fit1);
                end
            end
        else
            s_num=s_num+1;
        end  
        iter_num=iter_num+1;
    end   
    T=T*alpha;
    disp(T);
    if flag==0
        flag_num=flag_num+1;
    else
        flag_num=0;
    end
    disp(flag_num);
    flag=0;%若解的分数始终未上涨，标记该温度 
    toc
end
disp("end with simulation")
indmax=post_process(indmax,model,v,regnet,BPCYids,tar);
[objmax,~,~,target,growth,~,new_model]=OptRAM_obj(indmax,model,v,regnet,BPCYids,tar);
for z=1:numel(indmax.gene)
    fprintf(fid,'%d\t',indmax.gene(z));
end
for z=1:numel(indmax.gene)
    fprintf(fid,'%d\t',indmax.code(z));
end
fprintf(fid,'%f\t%f\t%f\n',objmax,target,growth);
fclose('all'); 
end

function individual=rand_ind(ngene,ntf,size)
    individual.gene=zeros(1,size);
    individual.code=zeros(1,size);
    %设置代谢基因
    mgenesize=round(rand*size);
    for i = 1:mgenesize
        individual.gene(i)=round(rand*(ngene-1))+1;
        individual.code(i)=ceil(rand*11)-6;%codenum可变
    end
    %设置调控基因
    tfsize=size-mgenesize;
    nmgene=ngene-ntf;
    for i = mgenesize+1:mgenesize+tfsize
        individual.gene(i)=round(rand*(ntf-1))+1;
        individual.code(i)=ceil(rand*11)-6;%codenum可变
        individual.gene(i)=individual.gene(i)+nmgene;
    end
end

function individual3=mutant(individual1,ngene,regnet)
    %突变种类，geneid改变，genecode改变,%%size改变
    size=length(individual1.code);
    mutant_id=ceil(rand*size);
    individual2.gene=individual1.gene;
    individual2.code=individual1.code;
    type=ceil(rand*5);
    ntf=numel(regnet.TFnames);
    nmgene=ngene-ntf;
    if size==1
        type=3;
    elseif size>=10%突变个体的基因数大于等于1，小于等于10
        type=4;
    end

    switch type
        case 1
            %换个基因
            individual2.gene(mutant_id)=round(rand*(ngene-1))+1;
        case 2
            %换个编码
            individual2.code(mutant_id)=ceil(rand*11)-6;
        case 3
            individual2.gene(size+1)=round(rand*(ngene-1))+1;%加个位点
            individual2.code(size+1)=ceil(rand*11)-6;
        case 4
            individual2.gene(mutant_id)=[];%减个位点
            individual2.code(mutant_id)=[];
        case 5%替换一个基因为调控因子
            if any(regnet.mgeneid==individual2.gene(mutant_id))
                x= regnet.mgeneid==individual2.gene(mutant_id);
                tf=find(regnet.mat(x,:)~=0);
                s=numel(tf);
                y=ceil(rand*s);
                if regnet.mat(x,tf(y)) > 0
                    individual2.gene(mutant_id)=tf(y)+nmgene;
                elseif  individual2.code(mutant_id)~=0
                    individual2.gene(mutant_id)=tf(y)+nmgene;
                    individual2.code(mutant_id)=-individual2.code(mutant_id);
                else
                    individual2.gene(mutant_id)=tf(y)+nmgene;
                    individual2.code(mutant_id)=5;
                end
            end
    end
    [individual3.gene,y]=unique(individual2.gene);
    individual3.code=individual2.code(y);
end

function individual=checkTF(individual1,regnet,n_mgene)
if ~isempty(individual1.gene(individual1.gene>n_mgene))
    del=[];
    TFs=individual1.gene(individual1.gene>n_mgene);
    TFs=TFs-n_mgene;
    %disp(TFs);
    a=regnet.mat(:,TFs);
    gene2del=find(sum(abs(a),2)>0);
    mgeneid=regnet.mgeneid;
    for i=1:numel(gene2del)
        x=find(individual1.gene==mgeneid(gene2del(i)));
        if ~isempty(x)
            del=[del x];
        end
    end
    individual=individual1;
    individual.code(del)=[];
    individual.gene(del)=[];
else
    individual=individual1;
end
end

