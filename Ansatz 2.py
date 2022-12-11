import math
import numpy as np

# Vorgaben
v_dd = 1.65
v_ss = 1.65     # -v_ss = -1.65V
i_q = 20e-6
u_offs = 10e-6
u_e_max = 1.55
u_e_min = -0.75
u_a_max = 0.55
u_a_min = -0.65
v_d0 = 2e4
f_t = 15e6
phi_r = 50
c_l = 1e-12

# Naturkonstanten
e = 1.6021918e-19
k_bolz = 1.3806226e-23
n_i = 1.45e16
epsilon_0 = 8.85419e-12
epsilon_r_si = 11.7
epsilon_r_sio2 = 3.9

# technologieabhängige Kennwerte
mu_0n = 0.0226  # Beweglichkeit n
mu_0p = 0.0120  # Beweglichkeit p
t_oxn = 7.06e-9  # Oxiddicke n
t_oxp = 6.97e-9  # Oxiddicke p
n_an = 3.476e22  # Majoritätsträgerkonzentration n
n_ap = 1.558e23  # Majoritätstragerkonzentration p
u_th0n = 0.5096  # Schwellspannung n
u_th0p = -0.6233  # Schwellspannung p
updelta_wn = 132.8e-9  # W-Unterätzung n
updelta_wp = 121.2e-9  # W-Unterätzung p
updelta_ln = 41.5e-9  # L-Unterätzung n
updelta_lp = 94.4e-9  # L-Unterätzung p

# berechnete Näherungen
gamma = 0.6667  # Substratgegensteuereffekt
phi = 0.8754  # Substratgegensteuereffekt
c_ox = 4.933e-3  # Flächen-Kapazitätsbelag
k_0n = 1.109e-4  # prozessspezifischer Steilheitsfaktor
k_0p = 0.592e-4  # prozessspezifischer Steilheitsfaktor
k_ln = 3.72e-14  # Kanallängenmodulationsfaktor
k_lp = 8.3e-15  # Kanallängenmodulationsfaktor
u_t = 25.86  # Temperaturspannung


def calc_uds_sat(i, s, ttype):
    if ttype == 'n':
        return math.sqrt(2*i/(k_0n*s))
    elif ttype == 'p':
        return math.sqrt(2*i/(k_0p*s))
    else:
        print('Error, no transistortype chosen (n/p)')
        return 0


def calc_uth(ubs, ttype):
    if ttype == 'n':
        return u_th0n + gamma * (math.sqrt(phi-ubs) - math.sqrt(phi))
    elif ttype == 'n':
        return u_th0p + gamma * (math.sqrt(phi-ubs) - math.sqrt(phi))
    else:
        print('Error, no transistortype chosen (n/p)')
        return 0


s1 = 3
f39 = 3

# Berechnung M10, M11, M12, M13
s13_temp = 0.5
ua = float(0)
while(ua <= 0.99 or ua >= 1.01):
    ugs_12 = calc_uds_sat(i=i_q, s=s13_temp, ttype='n') + u_th0n
    uth_13 = calc_uth(ubs=-ugs_12, ttype='n')
    uds_sat13 = ugs_12 - uth_13
    ua = ugs_12 + uds_sat13
    s13_temp += 0.02
s13 = s13_temp
print(s13)
print(ua)

# Berechnung M3, M9
s9_temp = 0.5
ua = float(0)
while(ua <= 1.09 or ua >= 1.11):
    ugs_9 = calc_uds_sat(i=i_q, s=s9_temp, ttype='p') + u_th0p
    uds_sat3 = calc_uds_sat(i=i_q, s=f39*s9_temp, ttype='p')
    ua = ugs_9 + uds_sat3
    s9_temp += 0.02
s9 = s9_temp
s3 = f39*s9_temp

# Berechnung M6
'''
s6_temp = 0.5
ue = float(0)
while(ue <=0.89 or ue >= 0.91):
    uds_sat6 = calc_uds_sat(i=i_q, s=s6_temp, ttype='n')
    uds_sat1 = calc_uds_sat(i=i_q, s=s1, ttype='n')
    ue = uds_sat6 + uds_sat1 + u_th0n
    s6_temp += 0.02
s6 = s6_temp
'''
s6 = 2 * i_q / (k_0n * (u_e_min + v_ss - calc_uds_sat(i=0.5*i_q, s=s1, ttype='n') - u_th0n)**2)

print('Mit M1 =', s1, 'ergibt sich: M3 =', s3, ', M6 =', s6, 'M9 = ', s9, ' und M13 =', s13)
