function [ pathScore ] = Score( graph )
%MACDEGREE 此处显示有关此函数的摘要
%   此处显示详细说明
%计算最大度
branchNum = length(graph.solution);
position = zeros(branchNum,1);
degree = zeros(graph.path_length,1);
for i = 1 : 1 : branchNum
    node = graph.branch(i,graph.branch_length(i)+1);
    for j = 1 : 1 : graph.path_length
        if (graph.path(j) == node)
            position(i) = j;
            degree(j) = degree(j) + 1;
            break;
        end;
    end;
end;
maxDegree = max(degree);

%计算改造位点在主路径上连接点之间的平均距离
sumDist = 0.0;
for i = 1 : 1 : branchNum
    for j = 1 : 1 : branchNum
        sumDist = sumDist + (position(i)-position(j));
    end;
end;
averageDist = (sumDist+1)/branchNum/(branchNum-1);

pathScore = averageDist*sum(graph.branch_length)/maxDegree;

end

