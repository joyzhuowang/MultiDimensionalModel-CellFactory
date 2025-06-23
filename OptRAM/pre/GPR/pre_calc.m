function [ new_model ] = pre_calc( model )
%PRE_CALC 此处显示有关此函数
%   此处显示详细说明
model = string_to_array(model);
postfix_expression = zeros(length(model.infix_long) , max(model.infix_long)+10);%infix_long表示rxn的grRule组成元素数量（括号，or，and，基因）
postfix_long = zeros(length(model.infix_long) , 1);
for i= 1 : 1 : length(model.infix_long)
    tmp = infix_to_postfix(model , i , 1 , model.infix_long(i));%删除括号，将基因在模型的位置标出来。实际上是对grRules的简化。
                                                                %例子：postfix_expression=[131,128,129,-1,-2,……],表示128和129两个基因关系为or，然后与131呈现and关系
                                                                      %postfix_expression=[97,98,-2,146,-1,……]，表示97和98关系为and，然后与146为or关系
                                                                      %postfix_expression=[350,13,-1,122,351,-1，-2……]，表示350和13是or关系，122和351是or关系，两组之间为and关系
    postfix_long(i) = tmp(length(tmp));%除去括号后，基因加上or和and的数量
    for j= 1 : 1 : postfix_long(i)
        postfix_expression(i,j) = tmp(j);
    end
    %disp(i);
end
model.postfix_long = postfix_long;
model.postfix_expression = postfix_expression;
new_model = model;
end

