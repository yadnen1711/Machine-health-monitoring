clear;
% SPECIFY DATA 
%%ѡ���ļ�����
stringname='100';
%ѡ����Ҫ������źŵ���
points=5000; 
PATH= '/Users/apple/Desktop/fyp/dataset'; % path, where data are saved
HEADERFILE= strcat(stringname,'.hea');      % header-file in text format
ATRFILE= strcat(stringname,'.atr');        % attributes-file in binary format
DATAFILE=strcat(stringname,'.dat');        % data-file
SAMPLES2READ=points;         % number of samples to be read
                      % in case of more than one signal:
                            % 2*SAMPLES2READ samples are read
   
% LOAD HEADER DATA --
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d',[1,3]);
nosig= A(1);  % number of signals
sfreq=A(2);   % sample rate of data
clear A;
for k=1:nosig
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);           % format; here only 212 is allowed
    gain(k)= A(2);              % number of integers per mV
    bitres(k)= A(3);            % bitresolution
    zerovalue(k)= A(4);         % integer value of ECG zero point
    firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
end
fclose(fid1);
clear A;

% LOAD BINARY DATA --
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
signald= fullfile(PATH, DATAFILE);            % data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);
M2H= bitshift(A(:,2), -4);
M1H= bitand(A(:,2), 15);
PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit
PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end
switch nosig
case 2
    M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
    M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
    TIME=(0:(SAMPLES2READ-1))/sfreq;
case 1
    M( : , 1)= (M( : , 1)- zerovalue(1));
    M( : , 2)= (M( : , 2)- zerovalue(1));
    M=M';
    M(1)=[];
    sM=size(M);
    sM=sM(2)+1;
    M(sM)=0;
    M=M';
    M=M/gain(1);
    TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
otherwise  % this case did not appear up to now!
    % here M has to be sorted!!!
    disp('Sorting algorithm for more than 2 signals not programmed yet!');
end
clear A M1H M2H PRR PRL;
fprintf(1,'\\n$> LOADING DATA FINISHED \n');

% LOAD ATTRIBUTES DATA -
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);
ATRTIME=[];
ANNOT=[];
sa=size(A);
saa=sa(1);
i=1;
while i<=saa
    annoth=bitshift(A(i,2),-2);
    if annoth==59
        ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
        ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
        i=i+3;
    elseif annoth==60
        % nothing to do!
    elseif annoth==61
        % nothing to do!
    elseif annoth==62
        % nothing to do!
    elseif annoth==63
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        i=i+hilfe/2;
    else
        ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        ANNOT=[ANNOT;bitshift(A(i,2),-2)];
    end
   i=i+1;
end
ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
ATRTIME(length(ATRTIME))=[];   % last line = EOF
clear A;
ATRTIME= (cumsum(ATRTIME))/sfreq;
ind= find(ATRTIME <= TIME(end));
ATRTIMED= ATRTIME(ind);
ANNOT=round(ANNOT);
ANNOTD= ANNOT(ind);

% DISPLAY DATA 
figure(1); clf, box on, hold on ;grid on ;
plot(TIME, M(:,1),'r');
if nosig==2
    plot(TIME, M(:,2),'b');
end
for k=1:length(ATRTIMED)
    text(ATRTIMED(k),0,num2str(ANNOTD(k)));
end
xlim([TIME(1), TIME(end)]);
xlabel('Time / s'); ylabel('Voltage / mV');
string=['ECG signal ',DATAFILE];
title(string);
fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');
% -
fprintf(1,'\\n$> ALL FINISHED \n');


M = M(1:5000,:);
points = 5000;



level=8; wavename='bior2.6';
%ecgdata=ECGsignalM1;
ecgdata=M;
figure(2);
plot(ecgdata(1:points));grid on ;axis tight;axis([1,points,-2,5]);
title('Original ECG signal');
%%%%%%%%%%����С���任8��
[C,L]=wavedec([ecgdata(:,1); ecgdata(:,2)],level,wavename);%modified with respect to 2015a version
%%%%%%%��ȡ�߶�ϵ����
A1=appcoef(C,L,wavename,1);
A2=appcoef(C,L,wavename,2);
A3=appcoef(C,L,wavename,3);
A4=appcoef(C,L,wavename,4);
A5=appcoef(C,L,wavename,5);
A6=appcoef(C,L,wavename,6);
A7=appcoef(C,L,wavename,7);
A8=appcoef(C,L,wavename,8);
%%%%%%%��ȡϸ��ϵ��
D1=detcoef(C,L,1);
D2=detcoef(C,L,2);
D3=detcoef(C,L,3);
D4=detcoef(C,L,4);
D5=detcoef(C,L,5);
D6=detcoef(C,L,6);
D7=detcoef(C,L,7);
D8=detcoef(C,L,8);
%%%%%%%%%%%%�ع�
A8=zeros(length(A8),1); %ȥ������Ư��,8���Ƶ��Ϣ
RA7=idwt(A8,D8,wavename);
RA6=idwt(RA7(1:length(D7)),D7,wavename);
RA5=idwt(RA6(1:length(D6)),D6,wavename);
RA4=idwt(RA5(1:length(D5)),D5,wavename);
RA3=idwt(RA4(1:length(D4)),D4,wavename);
RA2=idwt(RA3(1:length(D3)),D3,wavename);
D2=zeros(length(D2),1); %ȥ����Ƶ������2���Ƶ����
RA1=idwt(RA2(1:length(D2)),D2,wavename);
D1=zeros(length(D1),1);%ȥ����Ƶ������1���Ƶ����
DenoisingSignal=idwt(RA1,D1,wavename);
figure(3);
plot(DenoisingSignal);
title('Denoised ECG signal'); grid on; axis tight;axis([1,points,-2,5]);
ecgdata1 = ecgdata;
clear ecgdata;

%<pre name="code" class="cpp">
level=4;   
sr=360; 
%����ECG�ź�
%load ecgdata.mat;
%load ECGsignalM1.mat;
%load Rsignal.mat
mydata = DenoisingSignal;
ecgdata=mydata';
swa=zeros(4,points);%�洢��ò��Ϣ
swd=zeros(4,points);%�洢ϸ����Ϣ
signal=ecgdata(0*points+1:1*points); %ȡ���ź�

%��С��ϵ���ͳ߶�ϵ��
%��ͨ�˲��� 1/4 3/4 3/4 1/4
%��ͨ�˲��� -1/4 -3/4 3/4 1/4
%��������С��

for i=1:points-3
   swa(1,i+3)=1/4*signal(i+3-2^0*0)+3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
   swd(1,i+3)=-1/4*signal(i+3-2^0*0)-3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
end
j=2;
while j<=level
   for i=1:points-24
     swa(j,i+24)=1/4*swa(j-1,i+24-2^(j-1)*0)+3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
     swd(j,i+24)=-1/4*swa(j-1,i+24-2^(j-1)*0)-3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
   end
   j=j+1;
end
%����ԭ�źźͳ߶�ϵ����С��ϵ��
%figure(10);
%subplot(level+1,1,1);plot(ecgdata(1:points));grid on ;axis tight;
%title('ECG�ź���j=1,2,3,4�߶��µĳ߶�ϵ������');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swa(i,:));axis tight;grid on; xlabel('time');ylabel(strcat('a  ',num2str(i)));
%end
%figure(11);
%subplot(level,1,1); plot(ecgdata(1:points)); grid on;axis tight;
%title('ECG�źż�����j=1,2,3,4�߶��µĳ߶�ϵ����С��ϵ��');
%for i=1:level
%    subplot(level+1,2,2*(i)+1);
%    plot(swa(i,:)); axis tight;grid on;xlabel('time');
%    ylabel(strcat('a   ',num2str(i)));
%    subplot(level+1,2,2*(i)+2);
%    plot(swd(i,:)); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%����ԭͼ��С��ϵ��
%figure(12);
%subplot(level,1,1); plot(real(ecgdata(1:points)),'b'); grid on;axis tight;
%title('ECG�źż�����j=1,2,3,4�߶��µ�С��ϵ��');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swd(i,:),'b'); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%**************************************����������ֵ��**********************%
ddw=zeros(size(swd));
pddw=ddw;
nddw=ddw;
%С��ϵ���Ĵ���0�ĵ�
posw=swd.*(swd>0);
%б�ʴ���0
pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
%������ֵ��
pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);
%С��ϵ��С��0�ĵ�
negw=swd.*(swd<0);
ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
%������ֵ��
nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);
%������
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,points)=1;
%�����ֵ���ֵ,��������0
wpeak=ddw.*swd;
wpeak(:,1)=wpeak(:,1)+1e-10;
wpeak(:,points)=wpeak(:,points)+1e-10;

%�������߶��¼�ֵ��
%figure(13);
%for i=1:level
%    subplot(level,1,i);
%    plot(wpeak(i,:)); axis tight;grid on;
%ylabel(strcat('j=   ',num2str(i)));
%end
%subplot(4,1,1);
%title('ECG�ź���j=1,2,3,4�߶��µ�С��ϵ����ģ����ֵ��');

interva2=zeros(1,points);
intervaqs=zeros(1,points);
Mj1=wpeak(1,:);
Mj3=wpeak(3,:);
Mj4=wpeak(4,:);
%�����߶�3��ֵ��
figure(14);
plot (Mj3);
%title('�߶�3��С��ϵ����ģ����ֵ��');

posi=Mj3.*(Mj3>0);
%��������ֵ��ƽ��
thposi=(max(posi(1:round(points/4)))+max(posi(round(points/4):2*round(points/4)))+max(posi(2*round(points/4):3*round(points/4)))+max(posi(3*round(points/4):4*round(points/4))))/4;
posi=(posi>thposi/3);
nega=Mj3.*(Mj3<0);
%�󸺼���ֵ��ƽ��
thnega=(min(nega(1:round(points/4)))+min(nega(round(points/4):2*round(points/4)))+min(nega(2*round(points/4):3*round(points/4)))+min(nega(3*round(points/4):4*round(points/4))))/4;
nega=-1*(nega<thnega/4);
%�ҳ���0��
interva=posi+nega;
loca=find(interva);
for i=1:length(loca)-1
    if abs(loca(i)-loca(i+1))<80
       diff(i)=interva(loca(i))-interva(loca(i+1));
    else
       diff(i)=0;
    end
end
%�ҳ���ֵ��
loca2=find(diff==-2);
%������ֵ��
interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
%������ֵ��
interva2(loca(loca2(1:length(loca2))+1))=interva(loca(loca2(1:length(loca2))+1));
intervaqs(1:points-10)=interva2(11:points);
countR=zeros(1,1);
countQ=zeros(1,1);
countS=zeros(1,1);
mark1=0;
mark2=0;
mark3=0;
i=1;
j=1;
Rnum=0;
%*************************��������ֵ�Թ��㡣��R����ֵ������?y��QRS����㼰�յ�*******************%
while i<points
    if interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<points&interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%�󼫴�ֵ�ԵĹ����
       mark3= round((abs(Mj3(mark2))*mark1+mark2*abs(Mj3(mark1)))/(abs(Mj3(mark2))+abs(Mj3(mark1))));
%R������ֵ��
       R_result(j)=mark3-10;%Ϊ��-10������ֵ��
       countR(mark3-10)=1;
%���QRS�����
       kqs=mark3-10;
       markq=0;
     while (kqs>1)&&( markq< 3)
         if Mj1(kqs)~=0
            markq=markq+1;
         end
         kqs= kqs -1;
     end
  countQ(kqs)=-1;
  
%���QRS���յ�  
  kqs=mark3-10;
  marks=0;
  while (kqs<points)&&( marks<3)
      if Mj1(kqs)~=0
         marks=marks+1;
      end
      kqs= kqs+1;
  end
  countS(kqs)=-1;
  i=i+60;
  j=j+1;
  Rnum=Rnum+1;
 end
i=i+1;
end


%************************ɾ�����㣬����©���**************************%
num2=1;
while(num2~=0)
   num2=0;
%j=3,�����
   R=find(countR);
%�������
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%������R�����С��0.4RRmeanʱ,ȥ��ֵС��R��
for i=2:length(R)
    if (R(i)-R(i-1))<=0.4*RRmean
        num2=num2+1;
        if signal(R(i))>signal(R(i-1))
            countR(R(i-1))=0;
        else
            countR(R(i))=0;
        end
    end
end
end

num1=2;
while(num1>0)
   num1=num1-1;
   R=find(countR);
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%������R���������1.6RRmeanʱ,��С��ֵ,����һ�μ�?yR��
for i=2:length(R)
    if (R(i)-R(i-1))>1.6*RRmean
        Mjadjust=wpeak(4,R(i-1)+80:R(i)-80);
        points2=(R(i)-80)-(R(i-1)+80)+1;
%��������ֵ��
        adjustposi=Mjadjust.*(Mjadjust>0);
        adjustposi=(adjustposi>thposi/4);
%�󸺼���ֵ��
        adjustnega=Mjadjust.*(Mjadjust<0);
        adjustnega=-1*(adjustnega<thnega/5);
%������
        interva4=adjustposi+adjustnega;
%�ҳ���0��
        loca3=find(interva4);
        diff2=interva4(loca3(1:length(loca3)-1))-interva4(loca3(2:length(loca3)));
%�����м���ֵ��,�ҳ�����ֵ��
        loca4=find(diff2==-2);
        interva3=zeros(points2,1)';
        for j=1:length(loca4)
           interva3(loca3(loca4(j)))=interva4(loca3(loca4(j)));
           interva3(loca3(loca4(j)+1))=interva4(loca3(loca4(j)+1));
        end
        mark4=0;
        mark5=0;
        mark6=0;
    while j<points2
         if interva3(j)==-1;
            mark4=j;
            j=j+1;
            while(j<points2&interva3(j)==0)
                 j=j+1;
            end
            mark5=j;
%������
            mark6= round((abs(Mjadjust(mark5))*mark4+mark5*abs(Mjadjust(mark4)))/(abs(Mjadjust(mark5))+abs(Mjadjust(mark4))));
            countR(R(i-1)+80+mark6-10)=1;
            j=j+60;
         end
         j=j+1;
     end
    end
 end
end
%����ԭͼ�������?y���
%%%%%%%%%%%%%%%%%%%%%%%%%%?_ʼ��PT����
%��R����ǰ�Ĳ��üӴ����������СΪ100��Ȼ����㴰���ڼ���С�ľ���
%figure(20);
%plot(Mj4);
%title('j=4 ϸ��ϵ��'); hold on
%%%%%%%����ֱ����j=4ʱ��R������
Mj4posi=Mj4.*(Mj4>0);
%��������ֵ��ƽ��
Mj4thposi=(max(Mj4posi(1:round(points/4)))+max(Mj4posi(round(points/4):2*round(points/4)))+max(Mj4posi(2*round(points/4):3*round(points/4)))+max(Mj4posi(3*round(points/4):4*round(points/4))))/4;
Mj4posi=(Mj4posi>Mj4thposi/3);
Mj4nega=Mj4.*(Mj4<0);
%�󸺼���ֵ��ƽ��
Mj4thnega=(min(Mj4nega(1:round(points/4)))+min(Mj4nega(round(points/4):2*round(points/4)))+min(Mj4nega(2*round(points/4):3*round(points/4)))+min(Mj4nega(3*round(points/4):4*round(points/4))))/4;
Mj4nega=-1*(Mj4nega<Mj4thnega/4);
Mj4interval=Mj4posi+Mj4nega;
Mj4local=find(Mj4interval);
Mj4interva2=zeros(1,points);
for i=1:length(Mj4local)-1
    if abs(Mj4local(i)-Mj4local(i+1))<80
       Mj4diff(i)=Mj4interval(Mj4local(i))-Mj4interval(Mj4local(i+1));
    else
       Mj4diff(i)=0;
    end
end
%�ҳ���ֵ��
Mj4local2=find(Mj4diff==-2);
%������ֵ��
Mj4interva2(Mj4local(Mj4local2(1:length(Mj4local2))))=Mj4interval(Mj4local(Mj4local2(1:length(Mj4local2))));
%������ֵ��
Mj4interva2(Mj4local(Mj4local2(1:length(Mj4local2))+1))=Mj4interval(Mj4local(Mj4local2(1:length(Mj4local2))+1));
mark1=0;
mark2=0;
mark3=0;
Mj4countR=zeros(1,1);
Mj4countQ=zeros(1,1);
Mj4countS=zeros(1,1);
flag=0;
while i<points
    if Mj4interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<points&Mj4interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%�󼫴�ֵ�ԵĹ����,��R4�м�ֵ֮���������R�㡣

       mark3= round((abs(Mj4(mark2))*mark1+mark2*abs(Mj4(mark1)))/(abs(Mj4(mark2))+abs(Mj4(mark1))));
       Mj4countR(mark3)=1;
       Mj4countQ(mark1)=-1;
       Mj4countS(mark2)=-1;
       flag=1;
    end
    if flag==1
        i=i+200;
        flag=0;
    else
        i=i+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%�ҵ�MJ4��QRS�������ȱ�ٶ�R���©���?y�������?y���Ȳ�ȥϸ���ˡ�
%%%%%
%%%%%�Գ߶�4��R���?y�����ã���Ҫ�Ľ��ĵط�
%%%%%%
%figure(200);
%plot(Mj4);
%title('j=4');
%hold on;
%plot(Mj4countR,'r');
%plot(Mj4countQ,'g');
%plot(Mj4countS,'g');

%%%%%%%%%%%%%%%%%%%%%%%%%%Mj4������ҵ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rlocated=find(Mj4countR);
Qlocated=find(Mj4countQ);
Slocated=find(Mj4countS);
countPMj4=zeros(1,1);
countTMj4=zeros(1,1);
countP=zeros(1,1);
countPbegin=zeros(1,1);
countPend=zeros(1,1);
countT=zeros(1,1);
countTbegin = zeros(1,1);
countTend = zeros(1,1);
windowSize=100;
%%%%%%%%%%%%%%%%%%%%%%%P����?y%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rlocated Qlocated ���ڳ߶�4�µ�����
for i=2:length(Rlocated)
    flag=0;
    mark4=0;
    RRinteral=Rlocated(i)-Rlocated(i-1);
    for j=1:5:(RRinteral*2/3)
       % windowEnd=Rlocated(i)-30-j;
       windowEnd=Qlocated(i)-j;
        windowBegin=windowEnd-windowSize;
        if windowBegin<Rlocated(i-1)+RRinteral/3
            break;
        end
        %���ڵļ���Сֵ
        %windowposi=Mj4.*(Mj4>0);
        %windowthposi=(max(Mj4(windowBegin:windowBegin+windowSize/2))+max(Mj4(windowBegin+windowSize/2+1:windowEnd)))/2;
        [windowMax,maxindex]=max(Mj4(windowBegin:windowEnd));
        [windowMin,minindex]=min(Mj4(windowBegin:windowEnd));
        if minindex < maxindex &&((maxindex-minindex)<windowSize*2/3)&&windowMax>0.01&&windowMin<-0.1
            flag=1;
            mark4=round((maxindex+minindex)/2+windowBegin);
            countPMj4(mark4)=1;
            countP(mark4-20)=1;
            countPbegin(windowBegin+minindex-20)=-1;
            countPend(windowBegin+maxindex-20)=-1;
        end
        if flag==1
            break; 
        end
    end
    if mark4==0&&flag==0 %����û��P������R������1/3����ֵ-1
        mark4=round(Rlocated(i)-RRinteral/3);
        countP(mark4-20)=-1;
    end
end
 %plot(countPMj4,'g');       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%T����?y%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear windowBegin windowEnd maxindex minindex windowMax windowMin mark4 RRinteral;

windowSizeQ=100;
for i=1:length(Rlocated)-1
    flag=0;
    mark5=0;
    RRinteral=Rlocated(i+1)-Rlocated(i);
    for j=1:5:(RRinteral*2/3)
       % windowBegin=Rlocated(i)+30+j;
       windowBegin=Slocated(i)+j;
        windowEnd  =windowBegin+windowSizeQ;
        if windowEnd >Rlocated(i+1)-RRinteral/4
            break;
        end
        %%%%%�����ڵļ���Сֵ
        [windowMax,maxindex]=max(Mj4(windowBegin:windowEnd));
        [windowMin,minindex]=min(Mj4(windowBegin:windowEnd));
        if minindex < maxindex &&((maxindex-minindex)<windowSizeQ)&&windowMax>0.1&&windowMin<-0.1
            flag=1;
            mark5=round((maxindex+minindex)/2+windowBegin);
            countTMj4(mark5)=1;
            countT(mark5-20)=1;%�ҵ�T����ֵ��
            %%%%%ȷ��T����ʼ����յ�
            countTbegin(windowBegin+minindex-20)=-1;
            countTend(windowBegin+maxindex-20)=-1;
        end
        if flag==1
            break;
        end
    end
    if mark5==0 %����û��T������R���� ���1/3����ֵ-2
        mark5=round(Rlocated(i)+ RRinteral/3);
        countT(mark5)=-2;
    end
end
%plot(countTMj4,'g');
%hold off;        
figure(4);
plot(ecgdata(0*points+1:1*points)),grid on,axis tight,axis([1,points,-2,5]);
title('ECG signal with wave detection');
hold on
plot(countR,'r');
plot(countQ,'k');
plot(countS,'k');

for i=1:Rnum
    if R_result(i)==0
        break
    end
    plot(R_result(i),ecgdata(R_result(i)),'bo','MarkerSize',10,'MarkerEdgeColor','g');
end
plot(countP,'r');
plot(countT,'r');
plot(countPbegin,'k');
plot(countPend,'k');
plot(countTbegin,'k');
plot(countTend,'k');

hold off





locatedR = find(countR);
locatedQ = find(countQ);
locatedS = find(countS);
locatedP = find(countP);
locatedT = find(countT);
locatedPon = find(countPbegin);
locatedPoff = find(countPend);
keep locatedR locatedQ locatedS locatedP locatedT locatedPon locatedPoff ecgdata1;
ecgdata = ecgdata1;
clear ecgdata1;
length=size(locatedR,2);
[Q,R,S,countQ, countS] = matchQRS(locatedQ, locatedR, locatedS, 0);
[T,~] = matchRT(locatedT,R,190,0);
[P,~] = matchPR(locatedP, R, 200, 0);
[P,R,T];
%% 
P=ans(:,1);
T=ans(:,3);
%% 

[Pon,P,Poff] = matchPwave(P, locatedPon, locatedPoff, 80, 0);
location = [Pon P Poff Q R S T];


%ɾȥ�������
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
final_location=location(1:5000-count_empty,:);
ampR = zeros(5000-count_empty-1,1);
ampQ = zeros(5000-count_empty-1,1);
ampS = zeros(5000-count_empty-1,1);
ampP = zeros(5000-count_empty-1,1);
ampT = zeros(5000-count_empty-1,1);
QRSinterval = zeros(5000-count_empty-1,1);
RRinterval = zeros(5000-count_empty-1,1);
PRinterval = zeros(5000-count_empty-1,1);
QTinterval = zeros(5000-count_empty-1,1);
STinterval = zeros(5000-count_empty-1,1);
PPinterval = zeros(5000-count_empty-1,1);
Pduration = zeros(5000-count_empty-1,1);

for i = 1:4999-count_empty
    Pisempty = 0;
    Sisempty = 0;
    Qisempty = 0;
    Risempty = 0;
    Tisempty = 0;
    Ponisempty = 0;
    Poffisempty = 0;
    
    % peak value of PQRST
    
    if final_location(i,2) == 0
        ampP(i) = 0;
        Pisempty = 1;
    else
        ampP(i) = ecgdata(final_location(i,2));
    end
    
    if final_location(i,4) == 0
        ampQ(i) = 0;
        Qisempty = 1;
    else
        ampQ(i) = ecgdata(final_location(i,4));
    end
    
    if final_location(i,5) == 0
        ampR(i) = 0;
        Risempty = 1;
    else
        ampR(i) = ecgdata(final_location(i,5));
    end
    
    if final_location(i,6) == 0
        ampS(i) = 0;
        Sisempty = 1;
    else
        ampS(i) = ecgdata(final_location(i,6));
    end
    
    if final_location(i,7) == 0
        ampS(i) = 0;
        Tisempty = 1;
    else
        ampT(i) = ecgdata(final_location(i,7));
    end
    
    % PR interval
    
    if Pisempty || Risempty
        PRinterval(i) = 0;
    else
        PRinterval(i) = final_location(i,5)-final_location(i,2);
    end
    
    % QRS duration
    
    if Qisempty || Sisempty
        QRSinterval(i) = 0;
    else
        QRSinterval(i) = final_location(i,6)-final_location(i,4);
    end
    
    
    % RR interval
    
    if final_location(i+1,5)==0 || final_location(i,5)==0
        RRinterval(i) = 0;
    else
        RRinterval(i) = final_location(i+1,5)-final_location(i,5);
    end
    
    if RRinterval(i)>450 && RRinterval(i)<750
        RRinterval(i) = RRinterval(i)/2;
    elseif RRinterval(i)>750
        RRinterval(i) = RRinterval(i)/3;
    end
    
    % QT interval
    
    if Tisempty || Qisempty
        QTinterval(i) = 0;
    else
        QTinterval(i) = final_location(i,7)-final_location(i,5);
    end
    
    % ST interval
    
    if Sisempty || Tisempty
        STinterval(i) = 0;
    else
       STinterval(i) = final_location(i,7)-final_location(i,6); 
    end
    
    % P wave duration
    
    if Ponisempty || Poffisempty
        Pduration(i) = 0;
    else
        Pduration(i) = final_location(i,3)-final_location(i,1);
    end
    
    % PP interval
    
    if final_location(i+1,2)==0 || final_location(i,2)==0
        PPinterval(i) = 0;
    else
        PPinterval(i) = final_location(i+1,2)-final_location(i,2);
    end
    
    if PPinterval(i) >450 && PPinterval(i)<750
        PPinterval(i) = PPinterval(i)/2;
    elseif PPinterval(i)>750
        PPinterval(i) = PPinterval(i)/3;
    end
    
end
feature_set=[ampP,ampQ,ampR,ampS,ampT,QRSinterval,...
    RRinterval,PRinterval,QTinterval,STinterval,Pduration,PPinterval];
final_location(size(final_location,1),:)=[];
feature_set1=feature_set;
QRS_location1=final_location;
%% 

save('feature_set1','feature_set1');
save('QRS_location1','QRS_location1');
