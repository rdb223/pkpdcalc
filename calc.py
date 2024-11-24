from flask import Flask, request, jsonify
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import io
import base64

app = Flask(__name__)

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

# Endpoint to calculate PK/PD curve
@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.json
    drug = data['drug']
    dose = data['dose']
    frequency = data['frequency']
    organism = data['organism']
    weight = data['weight']
    renal_function = data['renal_function']

    # Connect to the database to retrieve breakpoint and drug properties
    conn = sqlite3.connect('./CLSI_breakpoints.db')
    cursor = conn.cursor()
    cursor.execute("SELECT mic, vd, half_life FROM drugs WHERE name = ?", (drug,))
    result = cursor.fetchone()
    conn.close()

    if not result:
        return jsonify({'error': 'Drug not found'}), 404

    mic, vd, half_life = result

    # Calculate elimination rate constant (ke)
    ke = np.log(2) / half_life

    # Define time points for a 24-hour period
    time_points = np.linspace(0, 24, 100)

    # Calculate concentrations over time
    concentrations = calculate_concentration(dose, frequency, vd, ke, time_points)

    # Plotting the concentration-time curve
    plt.figure()
    plt.plot(time_points, concentrations, label=f'{drug} Concentration')
    plt.axhline(y=mic, color='r', linestyle='--', label=f'MIC ({mic} μg/mL)')
    plt.xlabel('Time (hours)')
    plt.ylabel('Concentration (μg/mL)')
    plt.title(f'PK/PD Profile for {drug} ({dose} mg every {24//frequency}h)')
    plt.legend()
    plt.grid()

    # Save the plot to a bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    plot_base64 = base64.b64encode(buf.getvalue()).decode('utf8')
    buf.close()

    return jsonify({'concentrations': concentrations, 'plot': plot_base64})

if __name__ == '__main__':
    app.run(debug=True)
