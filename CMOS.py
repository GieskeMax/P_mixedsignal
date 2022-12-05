import math

# Vorgaben
v_dd = 1.65
v_ss = -1.65
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



def calc_i_d_active(beta, u_gs, u_th, u_ds):
    return beta * ((u_gs - u_th) * u_ds - 0.5 * u_ds ** 2)

def calc_i_d_sat(beta, u_gs, u_th, u_ds):
    return 0.5 * beta * (u_gs - u_th) ** 2

def calc_(beta, i_da):
    return math.sqrt(2 * beta * i_da)

def calc_g_mb(eta, g_m):
    return eta * g_m

def calc_g_ds(i_da, k_l, delta_l, l_eff):
    return i_da * k_l / (2 * delta_l * (l_eff - delta_l))

def calc_u_ds_sat(u_gs, u_th):
    return u_gs - u_th

def calc_beta_n(w_eff, l_delta_l):
    return mu_0n * c_ox * w_eff / l_delta_l

def calc_beta_p(w_eff, l_delta_l):
    return mu_0p * c_ox * w_eff / l_delta_l

def calc_c_ox(t_ox):
    return epsilon_0 * epsilon_r_sio2 / t_ox

def calc_u_th_n(u_bs):
    return u_th0n + gamma * (math.sqrt(phi - u_bs) - math.sqrt(phi))

def calc_u_th_p(u_bs):
    return u_th0p + gamma * (math.sqrt(phi - u_bs) - math.sqrt(phi))

def calc_v_d (g_m1, g_m2, g_mb2, r_ds1, r_ds2, r_d):
    return -(g_m1*r_ds1*(1+(g_m2+g_mb2)*r_ds2)*r_d)/(2*(r_ds2+r_ds1*(1+(g_m2+g_mb2)*r_ds2)+r_d))

def calc_r_cascode (r_ds, gm, gmb):
    return r_ds*(r_ds*(gm+gmb+1/r_ds)+1)



#Kaskode-Stromspiegel:
def calc_wl_s ():
    return 2*i_q/(k_0n*(u_a_min-v_ss-u_th0n)**2)

wl_10_12 = calc_wl_s()
print('M10, M12:', wl_10_12)

def calc_u_gs_10 (wl_s):
    return u_th0n+(math.sqrt(2*i_q/(k_0n*wl_s)))

u_gs_10 = calc_u_gs_10(wl_10_12)
print('U_gs10:',u_gs_10)

# Annahme: u_gs >= u_th -→ u_bs <= -u_th damit M10 in Sättigung ist
u_thk_11 = calc_u_th_n(-u_th0n)
print('u_thk_11:', u_thk_11)

def calc_wl_c (wl_s, u_gs, u_thk):
    return 2*i_q/(k_0n*(u_a_min-v_ss-u_gs-u_thk)**2)

wl_11_13 = calc_wl_c(wl_10_12, u_gs_10, u_thk_11)
print('M11, M13:', wl_11_13)



