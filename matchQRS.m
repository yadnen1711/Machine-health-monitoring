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
    %������������Ϊ�˷��������������״̬
    
    try
        tempR = inputR(countR);
    catch
        return
    end
    
    %QRS ��������ô���У�����6�����
    
    if tempS < tempQ && tempS < tempR
        %��S��Q��R��С��˵��S�ǵ���һ�㣬�����������S<Q<R��S<R<Q
        Q(i) = value_if_not_found;
        R(i) = value_if_not_found;
        S(i) = tempS;
        
        countS = countS + 1;
        
        haveQ(i) = 0;
        haveR(i) = 0;
        haveS(i) = 1;
        
    elseif tempQ < tempS && tempS < tempR
        %��Q<S<R, ֤���������û��R����Ҳ�п��������������S���������ܵ�ǰ�棬ϣ������������
        disp(strcat('Q ', string(countQ), '��', 'S ', string(countS),'֮��û��R'));
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
        %��R<Q<S, ֤���������û��Q��S����Ҳ�п��������������Q���������ܵ����棬ϣ������������
%         disp(strcat('R ', string(countR),'û��Q��S'));
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
        %R<S<Q,���ֱܷ���������������������
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
            
        
    elseif tempQ < tempR && tempR < tempS %������
        
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
        disp('�����ܳ��ֵ����')
    end
    
end



end