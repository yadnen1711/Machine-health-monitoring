function [Q,R,S,countQ, countR, countS, haveQ, haveR, haveS] = matchQRS(inputQ, inputR, inputS,value_if_not_found)

% if length(inputR)< length(inputQ) || length(inputR)< length(inputS)
%     disp('error, R should be the longest');
%     return
% end

% lR = length(inputR);
% lQ = length(inputQ);
% lS = length(inputS);

lgth = 5000;

haveQ = zeros(lgth,1);
haveR = zeros(lgth,1);
haveS = zeros(lgth,1);




R = zeros(lgth,1);
Q = zeros(lgth,1);
S = zeros(lgth,1);

countQ = 1;
countS = 1;
countR = 1;
for i = 1:lgth
    try
        tempQ = inputQ(countQ);
    catch
        tempQ = -inf;
    end
    
    try
        tempS = inputS(countS);
    catch
        tempS = inf;
    end
    %这俩做无限是为了方便进入最正常的状态
    
    try
        tempR = inputR(countR);
    catch
        return
    end
    
    %QRS 三个数怎么排列，共有6种情况
    
    if tempS < tempQ && tempS < tempR
        %若S比Q和R都小，说明S是单独一点，这种情况包含S<Q<R和S<R<Q
        Q(i) = value_if_not_found;
        R(i) = value_if_not_found;
        S(i) = tempS;
        
        countS = countS + 1;
        
        haveQ(i) = 0;
        haveR(i) = 0;
        haveS(i) = 1;
        
    elseif tempQ < tempS && tempS < tempR
        %若Q<S<R, 证明这次心跳没有R。（也有可能是这次心跳的S出了问题跑到前面，希望不是这样）
        disp(strcat('Q ', string(countQ), '和', 'S ', string(countS),'之间没有R'));
        disp(strcat('Q ', string(countQ), '=', string(tempQ)));
        disp(strcat('S ', string(countS), '=', string(tempS)));
        
        if isinf(tempQ)
            tempQ = value_if_not_found;
        end
        Q(i) = tempQ;
        S(i) = tempS;
        R(i) = value_if_not_found;
        
        countQ = countQ + 1;
        countS = countS + 1;
                
        haveQ(i) = 1;
        haveR(i) = 0;
        haveS(i) = 1;
        
    elseif tempR < tempQ && tempQ < tempS
        %若R<Q<S, 证明这次心跳没有Q和S。（也有可能是这次心跳的Q出了问题跑到后面，希望不是这样）
%         disp(strcat('R ', string(countR),'没有Q和S'));
%         disp(strcat('R ', string(countR), '=', string(tempR)));
%         disp(strcat('Q ', string(countQ), '=', string(tempQ)));
%         disp(strcat('S ', string(countS), '=', string(tempS)));
        
        Q(i) = value_if_not_found;
        S(i) = value_if_not_found;
        R(i) = tempR;
        
        countR = countR + 1;
        
        haveQ(i) = 0;
        haveR(i) = 1;
        haveS(i) = 0;
        
    elseif tempR < tempS && tempS < tempQ
        %R<S<Q,可能分别属于两个或者三个心跳
        if tempS <= tempR + 50
            Q(i) = value_if_not_found;
            S(i) = tempS;
            R(i) = tempR;
            countR = countR + 1;
            countS = countS + 1;
            
            haveQ(i) = 0;
            haveR(i) = 1;
            haveS(i) = 1;
        else
            Q(i) = value_if_not_found;
            S(i) = value_if_not_found;
            R(i) = tempR;
            countR = countR + 1;
            
            haveQ(i) = 0;
            haveR(i) = 1;
            haveS(i) = 0;
        end
            
        
    elseif tempQ < tempR && tempR < tempS %最正常
        
        if isinf(tempQ) && isinf(tempS)
            Q(i) = value_if_not_found;
            S(i) = value_if_not_found;
            R(i) = tempR;
            
            countR = countR + 1;
            
            haveQ(i) = 0;
            haveR(i) = 1;
            haveS(i) = 0;
            
        elseif isinf(tempQ)
            Q(i) = value_if_not_found;
            haveQ(i) = 0;
            if tempS <= tempR + 50
                R(i) = tempR;
                S(i) = tempS;
                
                countR = countR + 1;
                countS = countS + 1;
                
                haveR(i) = 1;
                haveS(i) = 1;
            else
                R(i) = tempR;
                S(i) = value_if_not_found;
                countR = countR + 1;
                
                haveR(i) = 1;
                haveS(i) = 0;
            end
            
        elseif isinf(tempS)
            S(i) = value_if_not_found;
            haveS(i) = 0;
            if tempQ >= tempR - 50
                Q(i) = tempQ;
                R(i) = tempR;
                countQ = countQ + 1;
                countR = countR + 1;
                
                haveQ(i) = 1;
                haveR(i) = 1;
            else
                Q(i) = temp(Q);
                R(i) = value_if_not_found;
                countQ = countQ + 1;
                haveQ(i) = 1;
                haveR(i) = 0;
            end
            
        else
            cnd1 = tempQ >= tempR - 50;
            cnd2 = tempS <= tempR + 50;
            
            if cnd1 && cnd2
                Q(i) = tempQ;
                R(i) = tempR;
                S(i) = tempS;
                countQ = countQ + 1;
                countR = countR + 1;
                countS = countS + 1;
                haveQ(i) = 1;
                haveR(i) = 1;
                haveS(i) = 1;
            elseif cnd1 && ~cnd2
                Q(i) = tempQ;
                R(i) = tempR;
                S(i) = value_if_not_found;
                countQ = countQ + 1;
                countR = countR + 1;
                haveQ(i) = 1;
                haveR(i) = 1;
                haveS(i) = 0;
            elseif ~cnd1 && cnd2
                Q(i) = tempQ;
                R(i) = value_if_not_found;
                S(i) = value_if_not_found;
                countQ = countQ + 1;
                haveQ(i) = 1;
                haveR(i) = 0;
                haveS(i) = 0;
            elseif ~cnd1 && ~cnd2
                Q(i) = tempQ;
                R(i) = value_if_not_found;
                S(i) = value_if_not_found;
                countQ = countQ + 1;
                haveQ(i) = 1;
                haveR(i) = 0;
                haveS(i) = 0;
            end
            
        end
        
        
    else
        disp('不可能出现的情况')
    end
    
end



end