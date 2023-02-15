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


# Step 1: Slew-rate
# keine Slew-rate als Parameter
i6 = i_q

# Step 2: Bias-Strom in Output-Kaskode
fi = 1.0            # Faktor von 1.2 bis 1.5
i2 = i3 = i8 = i9 = fi * i6    # vermeidet Null Strom in Kaskode bei Sperrung eines Eingangstransistors

# Step 3: Maximum Ausgangsspannung
usd2 = usd3 = usd8 = usd9 = 0.5 * (v_dd-u_a_max)

s8 = s9 = (2 * i8) / (k_0p * usd8**2)
s2 = s3 = (2 * i2) / (k_0p * usd2**2)
print(f"s2 = {s2}, s8 = {s8}")

# Step 4: Minimum Ausgangsspannung
uds10 = uds11 = uds12 = uds13 = 0.5 * (u_a_min + v_ss)
i10 = i11 = i12 = i13 = fi * i6
s10 = s11 = s12 = s13 = (2 * i10) / (k_0n * uds10**2)
print(f"s10 = {s10}")

# Step 5: Gain-Bandwith
# f_t = gm1 / c_l
s0 = s1 = ((2*np.pi*f_t)**2 * c_l**2) / (k_0n*i6)   # GB als Winkelfrequenz omega_t = 2*pi*f_t
print(f"s0 = {s0}")

# Step 6: Minimum Eingangsgleichtaktspannung
uds6 = u_e_min + v_ss - np.sqrt(i3 / (k_0n * s1)) - u_th0n
s6 = (2 * i6) / (k_0n * uds6**2)
print(f"s6 = {s6}")

s6_temp = 0.5
ue = float(0)
while(ue <= 0.89 or ue >= 0.91):
    uds_sat6 = calc_uds_sat(i=i6, s=s6_temp, ttype='n')
    uds_sat1 = calc_uds_sat(i=i6, s=s1, ttype='n')
    ue = uds_sat6 + uds_sat1 + u_th0n
    s6_temp += 0.02
print(f"s6_backgate = {s6_temp}")

# Step 7: Maximum Eingangsgleichtaktspannung
uds8_test = uds9_test = v_dd - u_e_max + u_th0n
s8_test = s9_test = (2 * i8) / (k_0p * uds8_test**2)
if (s8 >= s8_test):
    print(f"s8 checked, s8_test={s8_test} < s8")
elif (s8 < s8_test):
    print(f"s8 not checked, s8_test={s8_test}")

# Step 8: Differenzverstärkung
gm8 = gm9 = np.sqrt(2 * i8 * k_0p * s8)     # gm = sqrt( 2 * ida * beta)
# gds = ida / (uA + udsA + uds_sat)
#gds 8 = gds9 = i8


