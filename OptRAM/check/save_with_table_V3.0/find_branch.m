function [ branch, branch_length ] = find_branch( graph )

minScore = inf;
for i = 1 : 1 : graph.path_length
    [graph.branch, graph.branch_length] = SPFA(graph, i);
    tmp = Score(graph);
    if (tmp<minScore)
        minScore = tmp;
        branch = graph.branch;
        branch_length = graph.branch_length;
    end;
end;    
end