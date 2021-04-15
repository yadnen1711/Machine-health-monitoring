function [A,R,hasA,hasR, tot_count] = matchAR(annotation, reference, threshold, value_if_not_found)
%threshold 是定义"挨得近"


la = length(annotation);
lr = length(reference);

lgth = la + lr; %最长的情况

A = zeros(lgth,1);
R = zeros(lgth,1);

hasA = zeros(lgth,1);
hasR = zeros(lgth,1);

countA = 1;
countR = 1;
tot_count = lgth;
for i = 1:lgth
   try
       tempA = annotation(countA);
   catch
       tempA = nan;
   end
   
   try
       tempR = reference(countR);
   catch
       tempR = nan;
   end
   
   if isnan(tempA) && isnan(tempR)
       tot_count = i-1;
       break;
   elseif isnan(tempA) && ~isnan(tempR)
       R(i) = tempR;
       A(i) = value_if_not_found;
       
       countR = countR + 1;
       hasA(i) = 0;
       hasR(i) = 1;
       continue;
   elseif ~isnan(tempA) && isnan(tempR)
       R(i) = value_if_not_found;
       A(i) = tempA;
       
       countA = countA + 1;
       hasA(i) = 1;
       hasR(i) = 0;
       continue;
   end
   
   %都没有取完则进入正常状态
   %一共四种情况：是否挨得近 × 谁先谁后
   % nan与double相比较会返回 false
   
   is_close = abs(tempA - tempR) <= threshold;
   Alarger = tempA > tempR;
   
   if is_close 
       A(i) = tempA;
       R(i) = tempR;
       countR = countR + 1;
       countA = countA + 1;
       hasA(i) = 1;
       hasR(i) = 1;
   else
       if Alarger
           A(i) = value_if_not_found;
           R(i) = tempR;
           countR = countR + 1;
           hasA(i) = 0;
           hasR(i) = 1;
       else
           A(i) = tempA;
           R(i) = value_if_not_found;
           countA = countA + 1;
           hasA(i) = 1;
           hasR(i) = 0;
       end
   end
       
   
end

tot_count = tot_count - 1;

A = A(1:tot_count);
R = R(1:tot_count);
hasA = hasA(1:tot_count);
hasR = hasR(1:tot_count);



end