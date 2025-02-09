import numpy as np
import pandas as pd

# 材料和环境参数
E_eff = 2.016
E_maxwell = 2.016
E_kelvin = 80
eta1 = 1800
eta2 = 100
sigma_y = 0.31
cohension = 0.09
Jc_initial = 0.000050
Gf = 4
alpha = 0.1
num_points = 14
max_width = 0.1413  
tau_inital = 0.0003
height = 0.4
width_increase = 0.000596
depth_increase = 0.001775
Y = 1.12
rouh = 1000
rouh_soil = 1900
sigma_0 = 0.001
initial_width = 0.01
total_time = 500
dt = 1.0
num_steps = int(total_time / dt)
thickness = 0.01

#计算溃口的深度
def calculate_breach_depth (depth, depth_increase, dt, height) :
    return min(depth + depth_increase * dt, height)

# Burgers模型计算变形
def burgers_model(E_maxwell, E_kelvin, eta1, eta2,  dt, num_steps, tau_t, sigma_0):
    strain = np.zeros(num_steps)
    strain_dashpot1 = 0
    strain_dashpot2 = 0
    exp_factor = np.exp(-E_kelvin * dt / eta2)
    for i in range(1, num_steps):
        strain_spring1 = (tau_t + sigma_0) / E_maxwell
        strain_dashpot1 += ((tau_t + sigma_0) / eta1) * dt
        strain_spring2 = (tau_t + sigma_0) / E_kelvin
        strain_dashpot2 += (1 - exp_factor)
        strain[i] = strain_spring1 + strain_dashpot1 + strain_spring2 * strain_dashpot2
    return strain

# 计算应力强度因子K
def stress_intensity_factor(stress, current_crack_length):
    return stress * Y * np.sqrt(np.pi * current_crack_length)

# 计算J积分
def j_integral(K, E_eff):
    return (K ** 2) / E_eff

# 计算裂纹增量
def calculate_crack_increment(Jc_current, E_eff, stress, Y, cracks):
    return (Jc_current * E_eff) / (stress * Y ** 2 * np.pi )

# 计算应变能
def calculate_strain_energy(stress, strain, volume):
    return 0.5 * stress * strain * volume

# 保存结果到 Excel
def save_to_excel(time_steps, strains, crack_increments, crack_lengths, depths):
    output_data = np.column_stack((time_steps, strains, crack_increments, crack_lengths, depths))
    columns = ['时间 (s)'] + [f'应变值_点{i}' for i in range(num_points)] + ['裂纹增量 (m)', '裂纹长度 (m)', '深度 (m)']
    df = pd.DataFrame(output_data, columns=columns)
    df.to_excel('0.9m裂缝增长记录.xlsx', index=False)
    print("结果已保存到 '0.9m裂缝增长记录.xlsx'")

def main():
    widths = np.zeros(num_steps)
    strains = np.zeros((num_steps, num_points))
    crack_increments = np.zeros(num_steps)
    crack_lengths = np.zeros(num_steps)
    depths = np.zeros(num_steps)
    delta_Jp = 0
    time_steps = np.arange(num_steps) * dt
    current_crack_length = 0.0
    Jc_current = Jc_initial
    initial_crack_formed = False
    depth = 0.1
    sigma_w = 0
    sigma_u = 0
    cracks = 0
    breach_width = 0.1 + dt * width_increase * 0.01
    tau_t = tau_inital

    # 初始时计算中心区域的应变值
    center = num_points // 2
    half_initial_width = int(initial_width)
    center_start = max(center - half_initial_width, 0)
    center_end = min(center + half_initial_width, num_points - 1)
    initial_strain = burgers_model(E_maxwell, E_kelvin, eta1, eta2, dt, num_steps, tau_t, sigma_0)

    for p in range(center_start, center_end + 1):
        strains[0, p] = initial_strain[-1]

    for t in range(1, num_steps):
        current_half_width = half_initial_width + t
        center_start = max(center - current_half_width, 0)
        center_end = min(center + current_half_width, num_points - 1)
        depth = calculate_breach_depth(depth, depth_increase, dt, height)
        depths[t] = depth
        widths[t] = 0.1 + t * width_increase
        tau_t = tau_t * breach_width
        sigma_w = rouh * 9.81 * depth / 1000000
        sigma_u = rouh_soil * 9.81 * depth * 0.271 / 1000000
        stress = tau_t + sigma_w + sigma_u + sigma_0
        volume = thickness * depth * widths[t]

        if not initial_crack_formed:
            for p in range(num_points):
                if center_start <= p <= center_end:
                    strain = burgers_model(E_maxwell, E_kelvin, eta1, eta2, dt, num_steps, tau_t, sigma_0)
                    if t == 1:
                        strains[t, p] = strain[-1]
                    else:
                        strains[t, p] = strains[t - 1, p] + strain[-1]
                else:
                    if t > 1:
                        strains[t, p] = strains[t - 1, p]
         # 计算应变能
            strain_energy = calculate_strain_energy(stress, strains[t, center_start], volume)
            if strain_energy >=  Gf / 1000000 * depth * thickness and current_crack_length < depth:
                initial_crack_length = strain_energy / (2 * Gf / 1000000 * thickness)
                current_crack_length += initial_crack_length

                crack_increments[t] = initial_crack_length
                if current_crack_length > depth:
                    current_crack_length = depth
                initial_crack_formed = True
            K = stress_intensity_factor(stress, current_crack_length)
            J_actual = j_integral(K, E_eff)
            Jc_current = J_actual
        else:
            step_counter = t
            crack_increment = calculate_crack_increment(Jc_current, E_eff, stress, Y, cracks)
            current_crack_length =  current_crack_length + crack_increment
            cracks = (Jc_current * E_eff) / (stress * Y ** 2 * np.pi )
            crack_increments[step_counter] = crack_increment
            step_counter += 1
            K = stress_intensity_factor(stress, current_crack_length)

            # 计算新的临界J值
            delta_Jp = (K**2 / E_eff) * (1 + K**2 / (2 * current_crack_length * np.pi * (sigma_y ** 2)))
            Jc_updated = Jc_current + alpha * delta_Jp
            Jc_current = Jc_updated

            driving_moment = tau_t * depth  + (sigma_u + sigma_w) * max(0, depth / 3)
            resistance_moment = cohension*breach_width*thickness+2 * cohension * (depth - current_crack_length) * thickness
            if driving_moment > resistance_moment:
                print("水土压力力矩大于抵抗力矩，计算终止。")
                save_to_excel(time_steps[:t], strains[:t], crack_increments[:t], crack_lengths[:t], depths[:t])
                return
        crack_lengths[t] = current_crack_length

    # 准备输出数据
    output_data = np.column_stack((time_steps, strains, crack_increments[:num_steps], crack_lengths, depths))

    # 更新列名，包含深度变化
    columns = ['时间 (s)'] + [f'应变值_点{i}' for i in range(num_points)] + ['裂纹增量 (m)', '裂纹长度 (m)', '深度 (m)']

    # 将结果保存为 Excel 文件
    df = pd.DataFrame(output_data, columns=columns)
    df.to_excel('0.9m裂缝增长记录.xlsx', index=False)


if __name__ == "__main__":
    main()