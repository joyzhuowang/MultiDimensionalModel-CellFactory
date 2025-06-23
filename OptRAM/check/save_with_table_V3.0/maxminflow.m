function [ path_length, path, flow ] = maxminflow( graph )
%MAXMINFLOW �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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

