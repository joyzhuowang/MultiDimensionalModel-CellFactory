function [ path_length, path, flow ] = maxminflow( graph )
%MAXMINFLOW 此处显示有关此函数的摘要
%   此处显示详细说明
left = 0.0;
right = max(graph.worth);
while (right-left>0.00001)
    mid = (left+right)/2.0;
    [can_visit, ans_path, length] = check(graph, mid);
    if (can_visit==1)
        left = mid;
    else
        right = mid;
    end;
end;
[can_visit, ans_path, length] = check(graph, left);
path = ans_path;
path_length = length;
flow = left;
end

