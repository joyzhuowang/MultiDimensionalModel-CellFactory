function [ new_model ] = pre_calc( model )
%PRE_CALC �˴���ʾ�йش˺���
%   �˴���ʾ��ϸ˵��
model = string_to_array(model);
postfix_expression = zeros(length(model.infix_long) , max(model.infix_long)+10);%infix_long��ʾrxn��grRule���Ԫ�����������ţ�or��and������
postfix_long = zeros(length(model.infix_long) , 1);
for i= 1 : 1 : length(model.infix_long)
    tmp = infix_to_postfix(model , i , 1 , model.infix_long(i));%ɾ�����ţ���������ģ�͵�λ�ñ������ʵ�����Ƕ�grRules�ļ򻯡�
                                                                %���ӣ�postfix_expression=[131,128,129,-1,-2,����],��ʾ128��129���������ϵΪor��Ȼ����131����and��ϵ
                                                                      %postfix_expression=[97,98,-2,146,-1,����]����ʾ97��98��ϵΪand��Ȼ����146Ϊor��ϵ
                                                                      %postfix_expression=[350,13,-1,122,351,-1��-2����]����ʾ350��13��or��ϵ��122��351��or��ϵ������֮��Ϊand��ϵ
    postfix_long(i) = tmp(length(tmp));%��ȥ���ź󣬻������or��and������
    for j= 1 : 1 : postfix_long(i)
        postfix_expression(i,j) = tmp(j);
    end
    %disp(i);
end
model.postfix_long = postfix_long;
model.postfix_expression = postfix_expression;
new_model = model;
end

