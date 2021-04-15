inputTfunction [T,R] = matchRT(inputT, inputR, threshold, value_if_not_found)
%threshold 是定义"挨得近"


lt = size(inputT,2);
lr = size(inputR,1);

lgth = lt + lr; %最长的情况

T = zeros(lgth,1);
R = zeros(lgth,1);

hasT = zeros(lgth,1);
hasR = zeros(lgth,1);

countT = 1;
countR = 1;
for i = 1:lgth
   try
       tempT = inputT(countT);
   catch
       tempT = nan;
   end
   
   try
       tempR = inputR(countR);
   catch
       tempR = nan;
   end
   
   if isnan(tempT)
       break;
   end
   
   %都没有取完则进入正常状态
   %一共四种情况：是否挨得近 × 谁先谁后
   % nan与double相比较会返回 false
   
   is_close = abs(tempT - tempR) <= threshold;
   Tlarger = tempT > tempR;
   
   if is_close 
       T(i) = tempT;
       R(i) = tempR;
       countR = countR + 1;
       countT = countT + 1;
       hasT(i) = 1;
       hasR(i) = 1;
   else
       if Tlarger
           T(i) = value_if_not_found;
           R(i) = tempR;
           countR = countR + 1;
           hasT(i) = 0;
           hasR(i) = 1;
       else
           T(i) = tempT;
           R(i) = value_if_not_found;
           countT = countT + 1;
           hasT(i) = 1;
           hasR(i) = 0;
       end
   end
       
   
end
R(i:i+200) = inputR(countR:countR+200);

T = T(1:5000);
R = R(1:5000);


end