df = csvread('B:\Sem 5\EE321\14\Dataset_EEG.csv'); %#ok<CSVRD> 
X = df(1,:); %Input values

fs = 200;       %sampling frequency
Ts = 1/fs;

bw = (2*pi)*3;  %3Hz
K1 = 1;         %Passband Attenuation
K2 = 40;        %Stopband Attenuation
ripple = 1;     %Max possible ripple

w5 = (2*pi)*31; %Gamma Wave

w41 = (2*pi)*14;%Beta Wave
w42 = (2*pi)*30;

w31 = (2*pi)*9; %Alpha Wave
w32 = (2*pi)*14;

w21 = (2*pi)*4; %Theta Wave
w22 = (2*pi)*9;

w1 = (2*pi)*4;  %Delta Wave

order = orderFilter;
n5 = order.orderCheby_high(K2,ripple,w5-bw,w5);     %wr=w5-bw, Chebyshev Filter
n4 = order.orderButter_band(K2,K1,w41,w42,bw);      %Butterworth Filter
n3 = order.orderButter_band(K2,K1,w31,w32,bw);      %Butterworth Filter
n2 = order.orderButter_band(K2,K1,w21,w22,bw);      %Butterworth Filter
n1 = order.orderCheby_low(K2,ripple,w1+bw,w1);      %wr=w1+bw, Chebyshev Filter

chebyHP5 = chebyshevTF(n5, w5, 0, 0, ripple, "high");   %Cheby HP object
butterBP4 = butterTF(n4,0,w41,w42);                     %Butter BP object
butterBP3 = butterTF(n3,0,w31,w32);                     %Butter BP object
butterBP2 = butterTF(n2,0,w21,w22);                     %Butter BP object
chebyLP1 = chebyshevTF(n1, w1, 0, 0, ripple, "low");    %Cheby LP object

[a5,b5] = chebyHP5.cheby_TF();                          %TF creation of various objects
[a4,b4] = butterBP4.butter_band();
[a3,b3] = butterBP3.butter_band();
[a2,b2] = butterBP2.butter_band();
[a1,b1] = chebyLP1.cheby_TF();

a5
b5

blt = Bi_Linear_Transform(Ts);                          %creating object to calculate BLT
[A5,B5] = blt.calcBLT(a5,b5);                           %caps represent TF after BLT
[A4,B4] = blt.calcBLT(a4,b4);
[A3,B3] = blt.calcBLT(a3,b3);
[A2,B2] = blt.calcBLT(a2,b2);
[A1,B1] = blt.calcBLT(a1,b1);

Y5 = filter(A5,B5,X);                                   %filtering input using various filters
Y4 = filter(A4,B4,X);
Y3 = filter(A3,B3,X);
Y2 = filter(A2,B2,X);
Y1 = filter(A1,B1,X);

plt = plotGraph(Ts);                                    %object for plotting amplitude vs time
%plt.plotTime(Y5,"Filter 5: Gamma Wave");
%plt.plotTime(Y4,"Filter 4: Beta Wave");
%plt.plotTime(Y3,"Filter 3: Alpha Wave");
%plt.plotTime(Y2,"Filter 2: Theta Wave");
%plt.plotTime(Y1,"Filter 1: Delta Wave");
plt.plotFreq(Y2);