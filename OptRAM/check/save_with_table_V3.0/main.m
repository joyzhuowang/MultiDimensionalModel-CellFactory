function [ answer ] = main( model , start , goal , p , q ,  solution )
%MAIN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
graph = build(model , start , goal , p , q, solution);
[graph.path_length, graph.path, graph.path_flow] = maxminflow(graph);
[graph.branch, graph.branch_length] = find_branch(graph);
%����ͼ������ͼ������re_branch��branch��ȫһ�����Ͳ�����re_branch��
[graph.path_name, graph.branch_name] = change_to_name(graph);
graph.pathScore = Score(graph);
answer = graph;
end

