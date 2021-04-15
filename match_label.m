function [aR, rR] = match_label(annoR,refR,value_not_found)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
lgth = 5000;


haveannoR= zeros(lgth,1);
haverefR = zeros(lgth,1);


aR = zeros(lgth,1);
rR = zeros(lgth,1);

countaR = 1;
countrR = 1;
for i = 1:lgth
    try
        tempaR = annoR(countaR);
    catch
        tempaR = inf;
    end
    
    try
        temprR = refR(countrR);
    catch
        temprR = -inf;
    end
    %这俩做无限是为了方便进入最正常的状态
    
    if temprR < tempaR
       % r a...r a...r a
       % r a...r  ...r a
       % r a...r  ...r  ...r a
       % r a........
       % r a...  a.......
       % r a...r  ......
       if temprR-tempaR <= 25
           %normal
          aR(i) = tempaR;
          rR(i) = temprR;
          
          haveannoR(i) = 1;
          haverefR(i) = 1;
          
          countaR = countaR + 1;
          countrR = countrR + 1;
          
       elseif isinf(temprR) && isinf(tempaR)
          aR(i) = value_not_found;
          rR(i) = value_not_found;
          
          haveannoR(i) = 0;
          haverefR(i) = 0;
          
       elseif isinf(temprR) && ~isinf(tempaR)
          aR(i) = tempaR;
          rR(i) = value_not_found;
          
          haveannoR(i) = 1;
          haverefR(i) = 0;
          
          countaR = countaR + 1;
           
       elseif isinf(tempaR) && ~isinf(temprR)
          aR(i) = value_not_found;
          rR(i) = temprR;
          
          haveannoR(i) = 0;
          haverefR(i) = 1;
          
          countrR = countrR + 1;
          
       else
           %隔一个/两个
           aR(i) = value_not_found; 
           rR(i) = temprR;
           
           haveannoR(i) = 0;
           haverefR(i) = 1;
           
           countrR = countrR + 1;
       end
       
    end
       
    if tempaR <= temprR
        % a r...a r...a r
        % a r...a  ...a r
        % a r...a  ...a  ...a r
        if tempaR-temprR <= 25
           %normal
          aR(i) = tempaR;
          rR(i) = temprR;
          
          haveannoR(i) = 1;
          haverefR(i) = 1;
          
          countaR = countaR + 1;
          countrR = countrR + 1;
        else
           %隔一个/两个
           aR(i) = tempaR; 
           rR(i) = value_not_found;
           
           haveannoR(i) = 1;
           haverefR(i) = 0;
           
           countaR = countaR + 1;
       end
    end
    
    
end

rR(find(rR==-inf))=0;
aR(find(aR==inf))=0;

end

