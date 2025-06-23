function [ branch, branch_length ] = SPFA( graph, standard )
%FIND_BRANCH 此处显示有关此函数的摘要
%   此处显示详细说明
branch_length = zeros(length(graph.solution),1);
branch = zeros(length(graph.solution),graph.n);
for t = 1 : 1 : length(graph.solution)
    dist = zeros(1,graph.n);
    for i = 1 : 1 : graph.n
        dist(i) = graph.n*2;
    end;
    v = zeros(1,graph.n);
    q = zeros(1,graph.n+1);
    pre = zeros(1,graph.n);
    head = 0;
    tail = 1;
    q(1) = graph.solution(t);
    dist(q(1)) = 0;
    v(q(1)) = 1;
    while head<tail
        head = head + 1;
        now = q(head);
        v(now) = 0;
        st = graph.source_all_first(now);
        while (st>0)
            if (dist(now)+1<dist(graph.source_all_goal_node(st)))
                dist(graph.source_all_goal_node(st)) = dist(now)+1;
                pre(graph.source_all_goal_node(st)) = now;
                if (v(graph.source_all_goal_node(st))==0)
                    tail = tail + 1;
                    q(tail) = graph.source_all_goal_node(st);
                    v(graph.source_all_goal_node(st)) = 1;
                end;
            end;
            st = graph.source_all_next(st);
        end;
    end;
    mind = graph.n*2;
    standardd = graph.n*2;
    minn = -1;
    for i = 1 : 1 : graph.path_length
        if ((dist(graph.path(i))<mind) || ((dist(graph.path(i))==mind)&&(abs(i-standard)<standardd)))
            mind = dist(graph.path(i));
            standardd = abs(i-standard);
            minn = graph.path(i);
        end;
    end;
    if minn==-1
        branch(t,1) = -1;
        continue;
    end;
    branch_length(t)=mind;
    now = minn;
    tmp = dist(now)+1;
    while now~=0
        branch(t,tmp) = now;
        tmp = tmp - 1;
        now = pre(now);
    end;

end;
end

