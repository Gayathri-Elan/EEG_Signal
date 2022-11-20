%%%%    BUTTER

filt = butterTF(1,1*2*pi, 0 ,0);

[a,b] = filt.butter_low();
a
b


perform = Bi_Linear_Transform;
[num,den] = perform.calcBLT(a,b,0.005);

figure

Y = filter(num,den,row1);

df = csvread('B:\Sem 5\EE321\14\Dataset_EEG.csv');
row1 = df(1,:);
row2 = 0:0.005:(10-0.005);
% plot(row2,Y);
% xlabel('time');
% ylabel('amplitude');
% 
% hold on
% 
plt = plotGraph(1/200);
plt.plotFreq(Y);
plot(row2,Y);