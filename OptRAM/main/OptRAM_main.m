function [new_model,fitmax,indmax] =  OptRAM_main(model,v,regnet,BPCYids,loc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The main optimization algorithm with simulated annealing
% INPUTS
%model--��лģ��
%v--�ο���ֵ
%regnet--��������:matrix;mgene;TFnames;mgeneExp;TFexp
%BPCYids--BPCY�з�Ӧ��ID
%tar--Ŀ�꺯��ѡ�����
%
% OUTPUTS
% fitmax--���Ŀ�꺯��
% indmax--���Ŀ�꺯����Ӧ�������
%
tar=3;%tar��ʾtarget��ΪĿ��������С
%��ʼ��ģ���˻����
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
tseed=str2num(ctime((end-5):end));%�����������
rand('seed',tseed)%����ר��
ntf=numel(regnet.TFnames);
n_mgene=numel(model.genes);
ngene=ntf+n_mgene;
individual1=rand_ind(ngene,ntf,initial_size);%�������һ������ĺ���
individual11=checkTF(individual1,regnet,n_mgene);
fit1=OptRAM_obj(individual11,model,v,regnet,BPCYids,tar);
filename=[loc,'log.txt'];
fid=fopen(filename, 'w');
fitmax=0;
indmax.gene=0;
indmax.code=0;
flag_num=0;%����ķ���ʼ��δ���ǣ���Ǹ��¶�,����¼�������ٸ��¶�δ����
flag=0;
disp(ntf);
%while flag_num<100 || (fitmax==0 && flag_num<300)
while flag_num<15 || (fitmax==0 && flag_num<30)
    tic
    iter_num=1; 
    s_num=1;
    fprintf(fid,'T:%.5e\n',T);
    while iter_num<iter_max && s_num<s_max; 
        individual2=mutant(individual1,ngene,regnet);%����ͻ�亯��
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
        if exp(Deltafit/T)>R;%ͻ����������������ͻ��
            individual1=individual2; 
            fit1=fit2;
            x1=round(fit1*1000000)/1000000;%�ж����ֵ�ļ��㾫���趨ΪС�����6λ
            s_num=1;
            if x1>fitmax%��¼�������size��С�ĸ���
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
    flag=0;%����ķ���ʼ��δ���ǣ���Ǹ��¶� 
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
    %���ô�л����
    mgenesize=round(rand*size);
    for i = 1:mgenesize
        individual.gene(i)=round(rand*(ngene-1))+1;
        individual.code(i)=ceil(rand*11)-6;%codenum�ɱ�
    end
    %���õ��ػ���
    tfsize=size-mgenesize;
    nmgene=ngene-ntf;
    for i = mgenesize+1:mgenesize+tfsize
        individual.gene(i)=round(rand*(ntf-1))+1;
        individual.code(i)=ceil(rand*11)-6;%codenum�ɱ�
        individual.gene(i)=individual.gene(i)+nmgene;
    end
end

function individual3=mutant(individual1,ngene,regnet)
    %ͻ�����࣬geneid�ı䣬genecode�ı�,%%size�ı�
    size=length(individual1.code);
    mutant_id=ceil(rand*size);
    individual2.gene=individual1.gene;
    individual2.code=individual1.code;
    type=ceil(rand*5);
    ntf=numel(regnet.TFnames);
    nmgene=ngene-ntf;
    if size==1
        type=3;
    elseif size>=10%ͻ�����Ļ��������ڵ���1��С�ڵ���10
        type=4;
    end

    switch type
        case 1
            %��������
            individual2.gene(mutant_id)=round(rand*(ngene-1))+1;
        case 2
            %��������
            individual2.code(mutant_id)=ceil(rand*11)-6;
        case 3
            individual2.gene(size+1)=round(rand*(ngene-1))+1;%�Ӹ�λ��
            individual2.code(size+1)=ceil(rand*11)-6;
        case 4
            individual2.gene(mutant_id)=[];%����λ��
            individual2.code(mutant_id)=[];
        case 5%�滻һ������Ϊ��������
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

