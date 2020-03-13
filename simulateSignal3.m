clc;
clear;
t = reshape(1/48000:1/48000:10,480000,1);
d1 = 0.02;%音响和麦克风距离
d2 = 0.3;%d1是referrence signal传递的距离，d2是手移动后信号传递的距离
v = 0.015;
d3 = d2 + 2*v*t;
TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 17150;
receivedSignal = zeros(480000,1);
cosSignal = zeros(480000,8);
sinSignal = zeros(480000,8);
%8个频率下的signal
for freNum = 1:8
    fre = f + freNum * 350;
    phase1 = d1*2*pi*fre/SoundSpeed;
    phase2 = d3*2*pi*fre/SoundSpeed;
    receivedSignal = receivedSignal + reshape(awgn(cos(2*pi*fre*t+phase2),100) + awgn(cos(2*pi*fre*t+phase1),100),480000,1);
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000,1);
end
framelength = 1920;
DCvalueI = zeros(8,1);
DCvalueR = zeros(8,1);
sumy = 0;
sumxy = 0;
maxR = zeros(8,1);
maxI = zeros(8,1);
minR = zeros(8,1);
minI = zeros(8,1);
totDis = 0;
temp = zeros(120,1);
temp1 = zeros(120,1);
temp2 = zeros(120,1);
Hm = mfilt.cicdecim(16,17,3);
data1 = zeros(250,1); 
data2 = zeros(250,1);
t1 = reshape(0:119,120,1);
freCount = zeros(8,1);
numCount = 0;
var = zeros(8,1);
varsum = 0;
for i = 2:250
    framesignal=receivedSignal((i-1)*framelength+1:i*framelength);
    for freNum = 1:8
        cos1 = cosSignal(((i-1)*framelength+1:i*framelength),freNum);
        sin1 = sinSignal(((i-1)*framelength+1:i*framelength),freNum);
        InPhase = framesignal.*sin1;
        Quard = framesignal.*cos1;
        %将两路信号用CICfilter处理
        InPhaseFi = filter(Hm,InPhase);
        QuardFi = filter(Hm,Quard);

        InPhasePro = double(InPhaseFi);
        QuardPro = double(QuardFi);

        %将两路信号用LEVD算法处理并求出DC vector
        [ESCInPhase,DCvalueR(freNum),maxR(freNum),minR(freNum),freCount(freNum)] = LEVD2(InPhasePro,DCvalueR(freNum),freCount(freNum));
        [ESCQuard,DCvalueI(freNum),maxI(freNum),minI(freNum),freCount(freNum)] = LEVD2(QuardPro,DCvalueI(freNum),freCount(freNum));
    
        %除去信号中的static vector 并将两路信号转变成baseband 
        DVInPhase = InPhasePro - ESCInPhase;
        DVQuard = QuardPro - ESCQuard;
        DVBaseband = ReImToComp(DVInPhase, DVQuard);

        % 求出相位同时根据threshold继续对DC vector进行处理

        [ph,DCvalueR(freNum),DCvalueI(freNum)] = DCprocess(DVBaseband,maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum));
        
        fre = f + freNum * 350;
        %linear regression
        if freCount(freNum) == 1
            for a = 2:120
                temp(a,freNum) = ph(a) - ph(1); 
                temp(a,freNum) = temp(a,freNum)*SoundSpeed/(2*pi*fre);
            end
            sumy = sum(temp(:,freNum))+sumy;
            sumxy = sum(temp(:,freNum).*t1)+sumxy;
            numCount = numCount + 1;
        end
    end
    
    %find ignore frequency
    if numCount == 0
        data1(i) = 0;
       continue; 
    end
    deltax = 8*((120-1)*120*(2*120-1)/6-119*120*119/4);
    delta = (sumxy-sumy*(120-1)/2)/deltax*8/numCount;
    temp1 = delta*t;
    
    %get variance of each frequency
    for freNum = 1:8
        temp2 = temp1 - delta;
        var(freNum) = sqrt(sum(temp2));
        varsum = varsum + var(freNum);
    end
    
    for freNum = 1:8
        if var(freNum) > varsum/numCount
            freCount(freNum) = 0;    
        end
    end
    
    sumxy = 0;
    sumy = 0;
    numCount = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 1:8
       if freCount(freNum) == 1
           sumy = sum(temp(:,freNum))+sumy;
           sumxy = sum(temp(:,freNum).*t1)+sumxy;
           numCount = numCount + 1;
       end
    end
    if numCount == 0
        data1(i) = 0;
        continue;
    end

    delta = (sumxy-sumy*(120-1)/2)/deltax*8/numCount;
    distance = -delta*120/2;    
    data1(i) = distance; 
    totDis = totDis + 1.5*distance;
    data2(i) = totDis;
end
x = 1/25:1/25:10;
plot(x,-x*v,'b',x,data2,'r');
