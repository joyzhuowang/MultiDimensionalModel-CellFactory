function [ answer ] = calc_postfix( expression, long , position, mgene_p )
%CALC_POSTFIX 此处显示有关此函数的摘要
%   此处显示详细说明
stack = zeros(long(position) , 2);
top = 0;
answer = [0 1];
for i=1 : 1 : long(position)
    if expression(position,i)>=0 %如果expression大于等于0，则说明对应代谢基因，把代谢基因的mgene_p存入stack
        top = top+1;
        stack(top,1) = mgene_p(expression(position,i));%代表基因的mgene_p
        stack(top,2) = 1;%代表基因数量
    elseif expression(position,i) == -1
        top = top-2;
        tmp = ORc([stack(top+1,1) stack(top+1,2)],[stack(top+2,1) stack(top+2,2)]);%计算呈or关系的两个基因的mgene_p的算术平均
        top = top+1;
        stack(top,1) = tmp(1);
        stack(top,2) = tmp(2);
    elseif  expression(position,i) == -2
        top = top-2;
        tmp = ANDc([stack(top+1,1) stack(top+1,2)],[stack(top+2,1) stack(top+2,2)]);%计算结果为[min(stack(top+1,1)),stack(top+2,1)),1],也就是较小的mgene_p被保留
        top = top+1;
        stack(top,1) = tmp(1);
        stack(top,2) = tmp(2);
    end
    %disp(top);
    %disp(stack);
end
if top == 1
    answer = [stack(1,1) stack(1,2)];
end

