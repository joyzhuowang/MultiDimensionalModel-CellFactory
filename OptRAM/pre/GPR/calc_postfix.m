function [ answer ] = calc_postfix( expression, long , position, mgene_p )
%CALC_POSTFIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
stack = zeros(long(position) , 2);
top = 0;
answer = [0 1];
for i=1 : 1 : long(position)
    if expression(position,i)>=0 %���expression���ڵ���0����˵����Ӧ��л���򣬰Ѵ�л�����mgene_p����stack
        top = top+1;
        stack(top,1) = mgene_p(expression(position,i));%��������mgene_p
        stack(top,2) = 1;%�����������
    elseif expression(position,i) == -1
        top = top-2;
        tmp = ORc([stack(top+1,1) stack(top+1,2)],[stack(top+2,1) stack(top+2,2)]);%�����or��ϵ�����������mgene_p������ƽ��
        top = top+1;
        stack(top,1) = tmp(1);
        stack(top,2) = tmp(2);
    elseif  expression(position,i) == -2
        top = top-2;
        tmp = ANDc([stack(top+1,1) stack(top+1,2)],[stack(top+2,1) stack(top+2,2)]);%������Ϊ[min(stack(top+1,1)),stack(top+2,1)),1],Ҳ���ǽ�С��mgene_p������
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

