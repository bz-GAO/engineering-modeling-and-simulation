function cost = eval_cost(aim, volt , tempr , Q)
%aim代表要计算成本的目标，volt和tempr代表标定误差的标准数据
%Q用于计算成本，抽象传入方便一点
%本质是测对于aim这条链，测试集里的很多传感器的综合效果，所以每一个都要做一次
%---判断是否有价值---
%利用sum函数求出1的和也就是采取的点的个数
num_of_points = sum(aim);
if num_of_points < 4 || num_of_points > 10
    cost = 1e9;
    return;
end

%---找到基准点---

[Nsample , gene_len] = size(volt);%前者代表样本的个体数，因为这里为了运算成本低，一次读入大量数据
logic_index = logical(aim);%logical返回一些true和false组成的对应数组，也就是转换成逻辑形式，运行速度比find快一点
est_tempr =zeros(Nsample , gene_len);

%---插值---

for j = 1:Nsample

    pre_mea_volt = volt(j ,logic_index);
    pre_mea_tempr = tempr(j ,logic_index);
    est_tempr(j,:) = mycurvefitting(pre_mea_volt , pre_mea_tempr, volt(j,:));

end

%---计算误差和成本---

    errabs = abs(est_tempr - tempr);

    le0_4=(errabs<=0.4);
    le0_8=(errabs<=0.8);
    le1_2=(errabs<=1.2);
    le1_6=(errabs<=1.6);
    le2_0=(errabs<=2);
    g2_0=(errabs>2);

    sij=2*(le0_8-le0_4)+15*(le1_2-le0_8)+30*(le1_6-le1_2)+45*(le2_0-le1_6)+80000*g2_0;
    si = sum(sij, 2) + Q * num_of_points;%sum里的2代表按行求和
    cost = sum(si)/Nsample;
end
