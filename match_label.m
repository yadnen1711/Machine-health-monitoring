function [aR, rR] = match_label(annoR,refR,value_not_found)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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
    %������������Ϊ�˷��������������״̬
    
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
           %��һ��/����
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
           %��һ��/����
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

