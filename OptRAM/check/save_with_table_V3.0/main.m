function [ answer ] = main( model , start , goal , p , q ,  solution )
%MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
graph = build(model , start , goal , p , q, solution);
[graph.path_length, graph.path, graph.path_flow] = maxminflow(graph);
[graph.branch, graph.branch_length] = find_branch(graph);
%由于图是无向图，所以re_branch和branch完全一样，就不计算re_branch了
[graph.path_name, graph.branch_name] = change_to_name(graph);
graph.pathScore = Score(graph);
answer = graph;
end

