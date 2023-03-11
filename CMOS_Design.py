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
v_d0 = 2e4      # entspricht 86dB
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
        return np.sqrt(2*i/(k_0n*s))
    elif ttype == 'p':
        return -np.sqrt(2*i/(k_0p*s))
    else:
        print('Error, no transistortype chosen (n/p)')
        return 0


def calc_uth(ubs, ttype):
    if ttype == 'n':
        return u_th0n + gamma * (np.sqrt(phi-ubs) - np.sqrt(phi))
    elif ttype == 'p':
        return u_th0p - gamma * (np.sqrt(phi-ubs) - np.sqrt(phi))
    else:
        print('Error, no transistortype chosen (n/p)')
        return 0


def calc_wl(uds_sat, ttype):
    if ttype == 'n':
        return (2*i_q) / (k_0n*uds_sat**2)
    elif ttype == 'p':
        return (2*i_q) / (k_0p*uds_sat**2)
    else:
        print('Error, no transistortype chosen (n/p)')
        return 0


# --------------------------------------------------------------------------------------------------------------
# Dimensionierungsansatz nach CMOS analog circuit design (Allen, Hoilberg - 2012)

# Step 1: Slew-rate (entfällt, da keine Slew-rate als Parameter gegeben)

# Step 2: Bias-Strom in Output-Kaskode (entfällt da I_2 = I_3 = I_6)

# Step 3: Maximum output voltage
u_sd8 = 0.5 * (v_dd - u_a_max)
u_sd2 = 0.5 * (v_dd - u_a_max)
s8 = (2*i_q)/(k_0p*u_sd8**2)
s2 = (2*i_q)/(k_0p*u_sd2**2)

# Step 4: Minimum output voltage
# s10 = s11 = s12 = s13
u_ds10 = 0.5 * (u_a_min + v_ss)
s10 = (2*i_q)/(k_0n*u_ds10)

# Step 5: Gainbandwidth
# f_t = gm1 / c_l
s0 = ((2*np.pi*f_t)**2 * c_l**2)/(k_0n*i_q)

# Step 6: Minimum input CM-voltage
s6 = 2*i_q / (k_0n * (u_e_min + v_ss - np.sqrt(i_q/k_0n*s0) - u_th0n)**2)

# Step 7: Maximum input CM-voltage
# s8 muss >= als in Step 3 sein
s8_test = 2*i_q / (k_0p*(v_dd - u_e_max + u_th0n))

print(f"W/L Verhältnisse:")
print(f"Differenzpaar: {s0}; folded-cascode: {s2}; cascode-mirror: {s10}; "
      f"top_bank: {s8}; bot_bank: {s6}")

if(s8_test >= s8):
    print(f"Step 7 erfolgreich")
else:
    print(f"Step 7 nicht erfolgreich;\n s8 = {s8}; s8_test = {s8_test}")
print(f"\n")

# --------------------------------------------------------------------------------------------------------------
# Dimensionierungsansatz nach Praktikumsanleitung

# Annahme von W/L des Differenzpaars mit s0 = s1 = 20
s0 = 20

# Berechnung des W/L Verhältnis der unteren Strombank über die minimale CM-Eingangsspannung
uds_sat1 = calc_uds_sat(i_q, s0, ttype='n')
uds_sat6 = uds_sat6 = u_e_min + v_ss - uds_sat1 - u_th0n
s6 = calc_wl(uds_sat6, ttype='n')

# Berechnung der W/L Verhältnisse der oberen Strombank und gefalteten Kaskode über die maximale Ausgangsspannung
# Annahme s3 = 2 * s9
s9_min = 2 * i_q / (k_0p * (1.1)**2) * (1+1/np.sqrt(2))**2
s3_min = 3 * s9_min

# Berechnung des W/L Verhältnis des Kaskoden Stromspiegels über die minimale Eingangsspannung
uds_sat13 = (u_a_min + v_ss - u_th0n) / 2 + k_0n
s13 = calc_wl(uds_sat13, ttype='n')

# Test, ob das min. W/L Verhältnis der oberen Strombank der Vorgabe der maximalen Eingangsspannung standhält
u_gs9 = (calc_uds_sat(i_q, s9_min, ttype='p') - calc_uth(0, ttype='p'))
u_e_max_test = v_dd + calc_uth(uds_sat6, ttype='n') + u_gs9

# Berechnung des W/L Verhältnis des oberen Stromspiegels über die maximale Eingangsspannung
u_gs9 = -(v_dd - u_e_max + calc_uth(uds_sat6, ttype='n'))
uds_sat9 = u_gs9 - u_th0p
s9_max = calc_wl(uds_sat9, ttype='n')

# Test, ob das max. W/L Verhältnis der oberen Strombank (s9_2) der Vorgabe der maximalen Ausgangsspannung standhält
s3_max = 2 * s9_max
uds_sat3 = calc_uds_sat(i_q, s3_max, ttype='n')
u_a_max_test = v_dd - u_gs9 - uds_sat3


print(f"W/L Verhältnisse:")
print(f"Differenzpaar: {s0}; folded-cascode: {s3_min} bis {s3_max}; cascode-mirror: {s13}; "
      f"top_bank: {s9_min} bis {s9_max}; bot_bank: {s6}")

if(u_e_max_test >= u_e_max):
    print(f"Test über maximale Eingangspannung erfolgreich (s9)")
else:
    print(f"Test über maximale Eingangspannung nicht erfolgreich (s9)")

if(u_a_max_test >= u_a_max):
    print(f"Test über maximale Eingangspannung erfolgreich (s3)")
else:
    print(f"Test über maximale Eingangspannung nicht erfolgreich (s3)")


