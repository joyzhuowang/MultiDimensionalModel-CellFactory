function [ new_model ] = build( model , start , goal , p , q , solution)
%BUILD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   pΪflux
%   qΪ�޳���С����
n = length(model.mets);
m = length(model.rxns);
list = zeros(m,1);

for i = 1 : 1 : m
    new_model.worth(i) = abs(p(i));
    if (p(i)<0)
        new_model.rxns(i) = strcat('--',model.rxns(i));
    else
        new_model.rxns(i) = model.rxns(i);
    end;
end;

%���²���Ϊ����ԭ���罨���ڽӱ�
new_model.all_m = 0;

source_all_first = zeros(m,1);
source_all_next = zeros(100000,1);
source_all_goal_node = zeros(100000,1);

source_re_all_first = zeros(m,1);
source_re_all_next = zeros(100000,1);
source_re_all_goal_node = zeros(100000,1);

for i = 1 : 1 : n
    % ����Ҫ���Ƕ���ԭʼ�����Ƿ�ɾ��q���漰�ķ�Ӧ��
    if (q(i)==0) 
        continue;
    end;
    now = 0;
    for j = 1 : 1 : m
        %disp([i,j]);
        if model.S(i,j)~=0
            now = now + 1;
            list(now) = j;
        end;
    end;
    for j = 1 : 1 : m
        if model.S(i,j)~=0
            for k = 1 : 1 : now
                if (list(k)==j)
                    continue;
                end;
                
                new_model.all_m = new_model.all_m + 1;
                source_all_goal_node(new_model.all_m) = list(k);
                source_all_next(new_model.all_m) = source_all_first(j);
                source_all_first(j) = new_model.all_m;
                
                source_re_all_goal_node(new_model.all_m) = j;
                source_re_all_next(new_model.all_m) = source_re_all_first(list(k));
                source_re_all_first(list(k)) = new_model.all_m;
                
                new_model.all_m = new_model.all_m + 1;
                source_all_goal_node(new_model.all_m) = j;
                source_all_next(new_model.all_m) = source_all_first(list(k));
                source_all_first(list(k)) = new_model.all_m;
                
                source_re_all_goal_node(new_model.all_m) = list(k);
                source_re_all_next(new_model.all_m) = source_re_all_first(j);
                source_re_all_first(j) = new_model.all_m;
                %end;
            end;
        end;
    end;
    %disp(i);
end;

%���²���Ϊ����ʵ����ֵ�����ڽӱ�

for i = 1 : 1 : n
    for j = 1 : 1 : m
        model.S(i,j) = model.S(i,j) * p(j);
    end;
end;

first = zeros(m,1);
next = zeros(100000,1);
goal_node = zeros(100000,1);
new_model.m = 0;
re_first = zeros(m,1);
re_next = zeros(100000,1);
re_goal_node = zeros(100000,1);

for i = 1 : 1 : n
    if (q(i)==0)
        continue;
    end;
    now = 0;
    for j = 1 : 1 : m
        if model.S(i,j)<0
            now = now + 1;
            list(now) = j;
        end;
    end;
    for j = 1 : 1 : m
        if model.S(i,j)>0
            for k = 1 : 1 : now
                %d(j,list(k)) = 1;
                new_model.m = new_model.m + 1;
                goal_node(new_model.m) = list(k);
                next(new_model.m) = first(j);
                first(j) = new_model.m;
               
                re_goal_node(new_model.m) = j;
                re_next(new_model.m) = re_first(list(k));
                re_first(list(k)) = new_model.m;
            end;
        end;
    end;
    %disp(i);
end;

new_model.solution = solution;
new_model.n = m;

new_model.start = start;
new_model.goal = goal;

new_model.first = first;
new_model.next = next;
new_model.goal_node = goal_node;

new_model.re_first = re_first;
new_model.re_next = re_next;
new_model.re_goal_node = re_goal_node;

new_model.source_all_first = source_all_first;
new_model.source_all_next = source_all_next;
new_model.source_all_goal_node = source_all_goal_node;

new_model.source_re_all_first = source_re_all_first;
new_model.source_re_all_next = source_re_all_next;
new_model.source_re_all_goal_node = source_re_all_goal_node;
end

