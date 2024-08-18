function Y=complex_network(a)
load data1
load data2
load data3
load data4


T=[data1,data2,data3,data4];
x=[0 0 1 1];
Y=a;
net1 = complx_net_feat(minmax(T),[20 10 1],{'logsig','logsig','purelin'},'trainrp');
net1.trainParam.show = 1000;
net1.trainParam.lr = 0.04;
net1.trainParam.epochs = 7000;
net1.trainParam.goal = 1e-5;[net1] = train(net1,T,x);
save net1 net1
y = round(sim(net1,T));
end
