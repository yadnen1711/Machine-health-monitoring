locatedR=find(countR);
locatedQ=find(countQ);
locatedS=find(countS);

keep locatedR locatedQ locatedS ecgdata1;
ecgdata = ecgdata1;
clear ecgdata1;
length=size(locatedR,2);
[Q,R,S,countQ, countS] = matchQRS(locatedQ, locatedR, locatedS, 0);
QRS = [Q R S];
count_empty = 0;


for i=1:size(Q,1)
    if Q(i) == 0 && R(i) == 0
        temp1 = 1;
    else
        temp1 = 0;
    end
    if R(i) == 0 && S(i) == 0
        temp2 = 1;
    else
        temp2 = 0;
    end
    if temp1 && temp2
        count_empty = count_empty+1;
    end
end
QRS_location=QRS(1:5000-count_empty,:);
ampR = zeros(5000-count_empty-1,1);ampQ = zeros(5000-count_empty-1,1);ampS = zeros(5000-count_empty-1,1);
QRSinterval = zeros(5000-count_empty-1,1);RRinterval = zeros(5000-count_empty-1,1);


for i = 1:4999-count_empty
    Sisempty = 0;
    Qisempty = 0;
    Risempty = 0;
    if QRS_location(i,1) == 0
        ampQ(i) = 0;
        Qisempty = 1;
    else
        ampQ(i) = ecgdata(QRS_location(i,1));
    end
    
    if QRS_location(i,2) == 0
        ampR(i) = 0;
        Risempty = 1;
    else
        ampR(i) = ecgdata(QRS_location(i,2));
    end
    
    if QRS_location(i,3) == 0
        ampS(i) = 0;
        Sisempty = 1;
    else
        ampS(i) = ecgdata(QRS_location(i,3));
    end
    
    if Qisempty || Sisempty
        QRSinterval(i) = 0;
    else
        QRSinterval(i) = QRS_location(i,3)-QRS_location(i,1);
    end
    
    if QRS_location(i+1,2)==0 || QRS_location(i,2)==0
        RRinterval(i) = 0;
    else
        RRinterval(i) = QRS_location(i+1,2)-QRS_location(i,2);
    end
    
    if RRinterval(i)>450 && RRinterval(i)<750
        RRinterval(i) = RRinterval(i)/2;
        
    elseif RRinterval(i)>750
        RRinterval(i) = RRinterval(i)/3;
    end
end
feature_set=[ampR,ampQ,ampS,QRSinterval,RRinterval];
QRS_location(size(QRS_location,1),:)=[];
feature_set1=feature_set;
QRS_location1=QRS_location;
save('feature_set1');
save('QRS_location1');
