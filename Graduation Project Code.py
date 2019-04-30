PI = 3.1415926
p_p, t_p, v_p, i_p = 300000, 240, 0.7805, 2947.505  # 工作蒸汽 压力Pa，温度摄氏度，比容m3/kg，比焓 kJ/kg
p_h, t_h, v_h, i_h = 25000, 65, 6.2046, 2618.347  # 被引射蒸汽 压力Pa，温度 摄氏度，比容m3/kg，比焓 kJ/kg
p_c = 45000  # 压缩（混合）蒸汽 压力Pa
k_p, k_h = 1.3, 1.13  # 绝热指数  对于过热水蒸气k为1.3,干饱和水蒸气k为1.13
R_p, R_h = 463, 463  # 气体常数 J/kg-K

import math
a_pcr=math.sqrt(2*k_p/(k_p+1))*math.sqrt(p_p*v_p) #工作流体的临界速度
a_hcr=math.sqrt(2*k_h/(k_h+1))*math.sqrt(p_h*v_h) #引射流体的临界速度
theta_sqrt=a_hcr/a_pcr
theta=math.pow(theta_sqrt,2)
pi_ph = p_h / p_p #引射流体压力/工作流体压力，pi函数

def lamda(pi,k):
    a=math.pow(pi,(k-1)/k)
    b=(k+1)/(k-1)
    return math.sqrt((1-a)*b)

def q(lamda,k):
    a=(k+1)/2
    b=1/(k-1)
    c=(k-1)/(k+1)
    return math.pow(a,b)*lamda*math.pow(1-c*lamda*lamda,b)

def pi(lamda,k):
    a=(k-1)/(k+1)
    b=k/(k-1)
    return math.pow(1-a*lamda*lamda,b)

def lamda_q(q,k):
    a=(k-1)/(k+1)
    b=2*math.pow(q,k-1)/(k+1)
    c1=k-1
    c2=k+1
    x1,x2=0,1
    x=0.5*(x1+x2)
    y=math.pow(x,c1)-a*math.pow(x,c2)
    while abs(b-y)>0.0001:
        x = 0.5 * (x1 + x2)
        y = math.pow(x, c1) - a * math.pow(x, c2)
        if y<b:
            x1=x
        else:
            x2=x
    return x

lamda_ph=lamda(pi_ph,k_p)
q_ph=q(lamda_ph,k_p)
pi_hcr=pi(1,k_h)
pi_ps=pi_hcr*pi_ph
lamda_ps=lamda(pi_ps,k_p)
q_ps=q(lamda_ps,k_p)
null=p_p/p_c*q_ps

print("工作流体临界速度a_pcr为"+str(round(a_pcr,4))+"m/s")
print("引射流体临界速度a_hcr为"+str(round(a_hcr,4))+"m/s")
print("引射与工作流体临界速度的比值theta_sqrt为"+str(round(theta_sqrt,4)))
print("lamda_ph为"+str(round(lamda_ph,4)))
print("q_ph为"+str(round(q_ph,4)))
print("pi_hcr为"+str(round(pi_hcr,4)))
print("pi_ps为"+str(round(pi_ps,4)))
print("q_ps为"+str(round(q_ps,4)))
print("非工作区null为"+str(round(null,4)))
if null>1:
    print("因为p_p/p_c*q_ps大于1而q函数始终不大于1，所以q_c3为任何值时喷射器的工作都是可能的。"+"\n")

k1,k2=0.834,0.812
lamda_c3_array=[]
q_c3_array=[]
pi_c3_array=[]
q_h2_array=[]
lamda_h2_array=[]
pi_h2_array=[]
pi_c2_array=[]
u_array=[]

lamda_c3=1

for i in range(1,61,1):
    q_c3=q(lamda_c3,k_p)
    pi_c3=pi(lamda_c3,k_p)
    u_np2 = 1 / theta_sqrt * (p_h / p_c / q_c3 - p_h / p_p / q_ps) / (1 - p_h / p_c / q_c3)
    u=u_np2
    tmp=0.7
    u1=0.4
    while tmp-u1>0.001:
        q_h2 = u * theta_sqrt / (p_h / p_c * (1 + u * theta_sqrt) / q_c3 - p_h / p_p / q_ph)
        lamda_h2 = lamda_q(q_h2, k_h)
        pi_h2 = pi(lamda_h2, k_h)
        pi_c2 = p_h / p_c * pi_h2
        k3 = 1 + 0.9 * p_c / p_p * (pi_c3 - p_h / p_c) / (k_p * pi_hcr * lamda_c3 * q_ph)
        k4 = 1 + 0.9 * p_c / p_h * (pi_c3 - pi_c2) / (k_p * pi_hcr * lamda_c3 * q_h2)
        u1 = 1 / theta_sqrt * (k1 * lamda_ph - k3 * lamda_c3) / (k4 * lamda_c3 - k2 * lamda_h2)
        tmp=u
        if u1<tmp:
            u=u1

    lamda_c3_array.append(lamda_c3)
    q_c3_array.append(q_c3)
    pi_c3_array.append(pi_c3)
    q_h2_array.append(q_h2)
    lamda_h2_array.append(lamda_h2)
    pi_h2_array.append(pi_h2)
    pi_c2_array.append(pi_c2)
    u_array.append(u)

    lamda_c3-=0.01

n=1
for i in range(1,60,1):
    if u_array[n]<u_array[i]:
        n=i

u=round(u_array[n],4)
print("最大可达喷射系数为"+str(u)+"，对应的热力学函数值如下：")
lamda_c3=round(lamda_c3_array[n],4)
q_c3=round(q(lamda_c3,k_p),4)
pi_c3=round(pi(lamda_c3,k_p),4)
q_h2 = round(u * theta_sqrt / (p_h / p_c * (1 + u * theta_sqrt) / q_c3 - p_h / p_p / q_ph),4)
lamda_h2 = round(lamda_q(q_h2, k_h),4)
pi_h2 = round(pi(lamda_h2, k_h),4)
pi_c2 = round(p_h / p_c * pi_h2,4)
print("lamda_c3为："+str(lamda_c3))
print("q_c3为："+str(q_c3))
print("pi_c3为："+str(pi_c3))
print("q_h2为："+str(q_h2))
print("lamda_h2为："+str(lamda_h2))
print("pi_h2为："+str(pi_h2))
print("pi_c2为："+str(pi_c2)+"\n")

G_p=500/3600#工作气体空气流量kg/h
pi_pcr=pi_hcr
f_pcr=round(G_p*a_pcr/(k_p*pi_pcr*p_p),6)#计算喷嘴临界截面积
q_p1=q_ph
f_p1=round(f_pcr/q_p1,6)
w_p=13.8#管道中流速m/s,w_p=13.8时，管内径100mm
f_p=round(G_p*v_p/w_p,6)
f_3=round(f_pcr*p_p*(1+u)*theta_sqrt/(p_c*q_c3),6)
f_2=f_3
d_1=round(math.sqrt(4*f_p1/PI),4)#计算轴向尺寸
d_3=round(math.sqrt(4*f_3/PI),4)
if u>0.5:
    d_4 = round(1.55 * d_1 * (1 + u),4)
else:
    d_4 = round(3.4 * d_1 * math.sqrt(0.083 + 0.76 * u),4)
a=0.08#试验常数，弹性介质在0.07-0.09之间
if u>0.5:
    l_c1 = round((0.37 + u) / (4.4 * a) * d_1,4)
else:
    l_c1 = round(d_1 / 2 / a * (math.sqrt(0.083 + 0.76 * u) - 0.29),4)
l_c2=round((d_4-d_3)/2,4)
if d_3>d_4:
    l_c = l_c1
else:
    l_c = l_c1 + l_c2
l_k=round(8*d_3,4)#混合室（喉部）长度，混合室长度一般取6-10倍直径
w_c,rho_c=20,0.2308701 #w_c：出口速度m/s；rho_c：出口混合汽密度kg/m3，根据出口比焓和压力查表得
f_c=round(G_p*(1+u)/(rho_c*w_c),4)#出口面积
d_c=round(math.sqrt(4*f_c/PI),4)#出口直径
l_PI=round(6.5*(d_c-d_3),4)#扩散器长度，4-5度扩张角时，系数是6-7
print("几何尺寸如下：")
print("喷嘴临界截面积f_pcr为："+str(f_pcr)+"m2")
print("喷嘴出口截面工作流体面积f_p1为："+str(f_p1)+"m2")
print("喷嘴入口截面积f_p为："+str(f_p)+"m2")
print("混合室出口截面积f_3为："+str(f_3)+"m2")
print("混合室入口截面积f_2为："+str(f_2)+"m2")
d_p,d_pcr=round(math.sqrt(4*f_p/PI),4),round(math.sqrt(4*f_pcr/PI),4)
print("喷嘴入口直径d_p为："+str(d_p)+"m")
print("喷嘴临界直径d_pcr为："+str(d_pcr)+"m")
print("喷嘴出口直径d_1为："+str(d_1)+"m")
print("混合室出口直径d_3为："+str(d_3)+"m")
print("自由流束终截面直径d_4为："+str(d_4)+"m")
print("扩散器出口直径d_c为："+str(d_c)+"m")
print("喷嘴出口到混合室入口距离l_c为："+str(l_c)+"m")
print("混合室长度l_k为："+str(l_k)+"m")
print("扩散器长度l_PI为："+str(l_PI)+"m")