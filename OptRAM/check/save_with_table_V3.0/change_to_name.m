function [ name1 , name2 ] = change_to_name( graph )
%CHANGE_TO_NAME 此处显示有关此函数的摘要
%   此处显示详细说明
for i = 1 : 1 : graph.path_length
    name1(i) = graph.rxns(graph.path(i));
end;
for i = 1 : 1 : length(graph.solution)
    for j = 1 : 1 : graph.n
        if graph.branch(i,j)<=0
            break;
        end;
        name2(i,j) = graph.rxns(graph.branch(i,j));
    end;
end;
end

