import streamlit as st
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import io
import base64

# Function to calculate concentration over time (one-compartment model)
def calculate_concentration(dose, frequency, vd, ke, time_points):
    concentrations = []
    dosing_interval = 24 / frequency
    for t in time_points:
        concentration = 0
        for dose_time in np.arange(0, t, dosing_interval):
            concentration += (dose / vd) * np.exp(-ke * (t - dose_time)) if t >= dose_time else 0
        concentrations.append(concentration)
    return concentrations

# Streamlit interface
st.title("PK/PD Dosing App")

# User inputs
drug = st.text_input("Enter the drug name (e.g., cephalexin):")
dose = st.number_input("Enter the dose (mg):", min_value=0, step=100)
frequency = st.number_input("Enter the frequency (doses per day):", min_value=1, step=1)
weight = st.number_input("Enter patient weight (kg):", min_value=0.0, step=0.1)
renal_function = st.number_input("Enter renal function (e.g., CrCl in mL/min):", min_value=0.0, step=0.1)

if st.button("Calculate PK/PD Curve"):
    # Connect to the database to retrieve breakpoint and drug properties
    conn = sqlite3.connect('./CLSI_breakpoints.db')
    cursor = conn.cursor()
    cursor.execute("SELECT mic, vd, half_life FROM drugs WHERE name = ?", (drug,))
    result = cursor.fetchone()
    conn.close()

    if not result:
        st.error("Drug not found in the database.")
    else:
        mic, vd, half_life = result

        # Calculate elimination rate constant (ke)
        ke = np.log(2) / half_life

        # Define time points for a 24-hour period
        time_points = np.linspace(0, 24, 100)

        # Calculate concentrations over time
        concentrations = calculate_concentration(dose, frequency, vd, ke, time_points)

        # Plotting the concentration-time curve
        fig, ax = plt.subplots()
        ax.plot(time_points, concentrations, label=f'{drug} Concentration')
        ax.axhline(y=mic, color='r', linestyle='--', label=f'MIC ({mic} μg/mL)')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (μg/mL)')
        ax.set_title(f'PK/PD Profile for {drug} ({dose} mg every {24//frequency}h)')
        ax.legend()
        ax.grid()

        # Display the plot in Streamlit
        st.pyplot(fig)
