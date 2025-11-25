import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math

# Stała przenikalności magnetycznej
mu0 = 4 * math.pi * 1e-7

st.title("🔧 Symulator cewki powietrznej + magnesu neodymowego")
st.write("Aplikacja oblicza parametry cewki, pole magnetyczne i siłę działającą na magnes neodymowy.")

# ---------------------------------------------------------
# 1️⃣ PARAMETRY CEWKI
# ---------------------------------------------------------

st.header("1️⃣ Parametry cewki powietrznej")

D_outer = st.number_input("Średnica zewnętrzna cewki [m]", min_value=0.005, max_value=1.0, value=0.06)
D_inner = st.number_input("Średnica wewnętrzna cewki [m]", min_value=0.005, max_value=1.0, value=0.02)
coil_height = st.number_input("Wysokość cewki [m]", min_value=0.005, max_value=1.0, value=0.05)
wire_diameter = st.number_input("Średnica drutu (z izolacją) [m]", min_value=0.0002, max_value=0.01, value=0.001)
I = st.number_input("Prąd DC zasilania [A]", min_value=0.1, max_value=50.0, value=5.0)

# ---------------------------------------------------------
# 2️⃣ PARAMETRY MAGNESU
# ---------------------------------------------------------

st.header("2️⃣ Parametry magnesu neodymowego")

mag_diameter = st.number_input("Średnica magnesu [m]", min_value=0.005, max_value=0.1, value=0.018)
mag_height = st.number_input("Wysokość magnesu [m]", min_value=0.002, max_value=0.1, value=0.01)

# typowe wartości magnetyzacji magnesów neodymowych: 6e5 - 1.1e6 A/m
M = st.number_input("Magnetyzacja osiowa M [A/m]", min_value=1e5, max_value=2e6, value=8e5, step=1e5)

# ---------------------------------------------------------
# 3️⃣ PARAMETRY SYMULACJI
# ---------------------------------------------------------

st.header("3️⃣ Parametry symulacji")

samples = st.slider("Liczba próbek wzdłuż osi z (dokładność symulacji)", 200, 2000, 600)
z_margin = st.slider("Margines osi Z poza cewką [m]", 0.0, 0.1, 0.02)

# ---------------------------------------------------------
# OBLICZENIA CEWKI
# ---------------------------------------------------------

radial_thickness = (D_outer - D_inner) / 2
n_radial_layers = max(1, int(radial_thickness / wire_diameter))
axial_turns = max(1, int(coil_height / wire_diameter))

layer_radii = []
turns_per_layer = []

for layer in range(n_radial_layers):
    r = D_inner/2 + (layer + 0.5) * wire_diameter
    layer_radii.append(r)
    circumference = 2 * math.pi * r
    N_turns = max(1, int(circumference / wire_diameter))
    turns_per_layer.append(N_turns)

total_turns = sum(turns_per_layer) * axial_turns

# Gęstość prądu
A_wire = math.pi * (wire_diameter / 2)**2
J_current_density = I / A_wire

# Objętość cewki
coil_volume = math.pi * ((D_outer/2)**2 - (D_inner/2)**2) * coil_height
volumetric_AT = (I * total_turns) / coil_volume

# ---------------------------------------------------------
# MAGNES
# ---------------------------------------------------------

r_mag = mag_diameter / 2
V_mag = math.pi * r_mag**2 * mag_height
m_dipole = M * V_mag

# ---------------------------------------------------------
# GENEROWANIE GEOMETRII CEWKI (LISTA PĘTLI)
# ---------------------------------------------------------

loops = []
for i, r in enumerate(layer_radii):
    for k in range(axial_turns):
        z0 = -coil_height/2 + (k + 0.5) * (coil_height / axial_turns)
        loops.append((r, z0, turns_per_layer[i]))

# ---------------------------------------------------------
# FUNKCJE PÓL MAGNETYCZNYCH
# ---------------------------------------------------------

def Bz_loop(z, a, z0, I_eff):
    zz = z - z0
    return mu0 * I_eff * a*a / (2 * (a*a + zz*zz)**1.5)

z_vals = np.linspace(-coil_height/2 - z_margin,
                     coil_height/2 + z_margin,
                     samples)

Bz_total = np.zeros_like(z_vals)

for (a, z0, Nturns) in loops:
    Bz_total += Bz_loop(z_vals, a, z0, I * Nturns)

# Pochodna ∂B/∂z
dBz_dz = np.gradient(Bz_total, z_vals)

# ---------------------------------------------------------
# SIŁA NA MAGNES
# ---------------------------------------------------------

def force_on_magnet(zc, points=20):
    zs = np.linspace(zc - mag_height/2, zc + mag_height/2, points)
    dz = zs[1] - zs[0]
    m_prime = M * math.pi * r_mag**2
    return np.sum(m_prime * np.interp(zs, z_vals, dBz_dz) * dz)

F_vals = np.array([force_on_magnet(z) for z in z_vals])

# ---------------------------------------------------------
# WYNIKI
# ---------------------------------------------------------

st.header("4️⃣ Wyniki obliczeń")

st.write(f"**Liczba warstw radialnych:** {n_radial_layers}")
st.write(f"**Przybliżona liczba zwojów:** {total_turns}")
st.write(f"**Gęstość prądu J:** {J_current_density:.2e} A/m²")
st.write(f"**Ampere-turns (I*N):** {(I * total_turns):.2e} A")
st.write(f"**Gęstość amperozwojowa (I*N)/V:** {volumetric_AT:.2e} A/m³")
st.write(f"**Moment dipolowy magnesu m:** {m_dipole:.2e} A·m²")

st.write("---")
st.write(f"**Maksymalna siła:** {np.max(F_vals):.5f} N")
st.write(f"**Minimalna siła:** {np.min(F_vals):.5f} N")

# ---------------------------------------------------------
# WYKRESY
# ---------------------------------------------------------

st.header("5️⃣ Wykresy pola magnetycznego i siły")

fig1, ax1 = plt.subplots()
ax1.plot(z_vals, Bz_total)
ax1.set_title("Pole magnetyczne B(z) na osi cewki")
ax1.set_xlabel("z [m]")
ax1.set_ylabel("Bz [T]")
ax1.grid(True)
st.pyplot(fig1)

fig2, ax2 = plt.subplots()
ax2.plot(z_vals, F_vals)
ax2.set_title("Siła F(z) działająca na magnes")
ax2.set_xlabel("z [m]")
ax2.set_ylabel("F [N]")
ax2.grid(True)
st.pyplot(fig2)

st.success("Symulacja zakończona!")
