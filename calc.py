import streamlit as st
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import io
import base64
import requests
import xml.etree.ElementTree as ET

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

# Function to query a drug database for pharmacokinetic information
def get_pharmacokinetics(drug_name):
    # Example using the PubChem API to get drug information
    api_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/xml'
    response = requests.get(api_url)

    if response.status_code == 200:
        root = ET.fromstring(response.content)
        half_life = None
        vd = None

        for info in root.iter('Information'):  # Searching for half-life information
            for comment in info.iter('String'): 
                if 'half-life' in comment.text.lower():
                    half_life = float(comment.text.split()[0])  # Assuming the value is the first part of the text
                    break

        if half_life:
            vd = 0.7  # Placeholder value for volume of distribution, update this if needed

        return half_life, vd
    else:
        return None, None

# Streamlit interface
st.title("PK/PD Dosing App")

# User inputs
drug = st.text_input("Enter the drug name (e.g., cephalexin):")
dose = st.number_input("Enter the dose (mg):", min_value=0, step=100)
frequency = st.number_input("Enter the frequency (doses per day):", min_value=1, step=1)
weight = st.number_input("Enter patient weight (kg):", min_value=0.0, step=0.1)
renal_function = st.number_input("Enter renal function (e.g., CrCl in mL/min):", min_value=0.0, step=0.1)

if st.button("Calculate PK/PD Curve"):
    # Query pharmacokinetic data
    half_life, vd = get_pharmacokinetics(drug)

    if not half_life or not vd:
        st.error("Drug information not found in the external database.")
    else:
        # Calculate elimination rate constant (ke)
        ke = np.log(2) / half_life

        # Define time points for a 24-hour period
        time_points = np.linspace(0, 24, 100)

        # Calculate concentrations over time
        concentrations = calculate_concentration(dose, frequency, vd, ke, time_points)

        # Plotting the concentration-time curve
        fig, ax = plt.subplots()
        ax.plot(time_points, concentrations, label=f'{drug} Concentration')
        ax.axhline(y=1, color='r', linestyle='--', label=f'MIC (placeholder)')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (μg/mL)')
        ax.set_title(f'PK/PD Profile for {drug} ({dose} mg every {24//frequency}h)')
        ax.legend()
        ax.grid()

        # Display the plot in Streamlit
        st.pyplot(fig)
