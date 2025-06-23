function [ can_visit, ans_path, length ] = check( graph, test_worth )
%CHECK 此处显示有关此函数的摘要
%   此处显示详细说明
ans_path = zeros(graph.n, 1);
length = 0;
if (graph.worth(graph.start)<test_worth)
    can_visit = 0;
    return;
end;
v = zeros(graph.n, 1);
q = zeros(graph.n ,1);
from = zeros(graph.n, 1);
head = 0;
tail = 1;
q(head+1) = graph.start;
v(graph.start) = 1;
while (head<tail)
    head = head + 1;
    now = q(head);
    t = graph.first(now);
    while (t>0)
        if ((v(graph.goal_node(t))==0)&&(graph.worth(graph.goal_node(t))>=test_worth))
            tail = tail + 1;
            q(tail) = graph.goal_node(t);
            from(q(tail)) = now;
            v(q(tail)) = 1;
        end;
        t = graph.next(t);
    end;
end;
if (v(graph.goal)==0)
    can_visit = 0;
    return;
end;
can_visit = 1;
len = 0;
now = graph.goal;
while (now~= 0)
    len = len + 1;
    now = from(now);
end;
length = len;
now = graph.goal;
while (now~= 0 )
    ans_path(len) = now;
    len = len -1;
    now = from(now);
end;
end

