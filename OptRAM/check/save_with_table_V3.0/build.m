function [ new_model ] = build( model , start , goal , p , q , solution)
%BUILD 此处显示有关此函数的摘要
%   此处显示详细说明
%   p为flux
%   q为剔除的小分子
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

%以下部分为对于原网络建立邻接表
new_model.all_m = 0;

source_all_first = zeros(m,1);
source_all_next = zeros(100000,1);
source_all_goal_node = zeros(100000,1);

source_re_all_first = zeros(m,1);
source_re_all_next = zeros(100000,1);
source_re_all_goal_node = zeros(100000,1);

for i = 1 : 1 : n
    % 这里要考虑对于原始网络是否删除q中涉及的反应物
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

%以下部分为考虑实际流值后建立邻接表

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

