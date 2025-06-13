import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import matplotlib.patches as mpatches

# Global parameters
L = 1.0
C = 0.00001
R1 = 1000.0
R2 = 1000.0
T = 0.2
dT = 0.001
amplitude = 1.0
frequency = 10
signal_type = 'step' # default
duty_cycle = 50.0  # %
harmonic_amplitudes = [0.3, 0.6, 0.9] # V
harmonic_frequencies = [10, 30, 50] # Hz
harmonic_phases = [0, 0, 0] # rad
input_old = [0]
output_old = [0]
time_old = [0]
param_old = [L, C, R1, R2]



def generate_signal():
    global L, C, R1, R2, T, dT, amplitude, frequency, signal_type, duty_cycle, harmonic_amplitudes, harmonic_frequencies, harmonic_phases
    
    length = int(T / dT)
    U = []
    
    if signal_type == "pulse":
        U = [amplitude]
        for i in range(length-1):
            U.append(0)
    elif signal_type == "step":
        for i in range(length):
            U.append(amplitude)
    elif signal_type == "sine":
        for i in range(length):
            t = i * dT
            U.append(amplitude * np.sin(2 * np.pi * frequency * t))
    elif signal_type == "square":
        for i in range(length):
            t = i * dT
            period = 1 / frequency
            time_in_period = (t % period) / period
            duty_ratio = duty_cycle / 100.0  # Konwersja z procentów na ułamek
            if time_in_period <= duty_ratio:
                U.append(amplitude)
            else:
                U.append(0)
    elif signal_type == "saw":
        for i in range(length):
            t = i * dT
            period = 1 / frequency
            time_in_period = (t % period) / period
            U.append(amplitude * (2 * time_in_period - 1))
    elif signal_type == "harmonic":
        for i in range(length):
            t = i * dT
            signal_value = 0
            for j in range(len(harmonic_amplitudes)):
                if j < len(harmonic_frequencies) and j < len(harmonic_phases):
                    signal_value += (harmonic_amplitudes[j] * 
                                   np.sin(2 * np.pi * harmonic_frequencies[j] * t + harmonic_phases[j]))
            U.append(amplitude * signal_value)
    
    return U

def calc_coefficients():
    global L, C, R1, R2, dT
    
    a = -dT / L
    b = dT / C
    c = -dT / (R1 * C) - dT / (R2 * C)
    d = dT / L
    
    g1 = L * C
    g2 = (L * (R1 + R2)**2 - 2 * (R1**2) * (R2**2) * C) / \
         ((R1**2) * (R2**2) * (C**2) * L)
    g3 = 1 / ((L**2) * (C**2))
    
    p1 = 1 / (R1 * C) + 1 / (R2 * C)
    p2 = 1 / (L * C)
    
    resonant_freq = np.sqrt(1 / (L * C))
    
    return [a, b, c, d, g1, g2, g3, p1, p2, resonant_freq]

K_old = calc_coefficients()

def runge_kutta(U): # Runge-Kutta second-order method
    K = calc_coefficients()
    length = len(U)
    
    x1 = [0]  # Current through inductor
    x2 = [0]  # Voltage across capacitor
    y = []    # Output voltage
    time = []
    
    for i in range(length - 1):
        x1_curr = x1[i]
        x2_curr = x2[i]
        u_curr = U[i]
        u_next = U[i + 1]
        
        # k1 calculations
        k1_x1 = K[0] * x2_curr + K[3] * u_curr
        k1_x2 = K[1] * x1_curr + K[2] * x2_curr
        
        xp1 = x1_curr + k1_x1
        xp2 = x2_curr + k1_x2
        
        # k2 calculations
        k2_x1 = K[0] * xp2 + K[3] * u_next
        k2_x2 = K[1] * xp1 + K[2] * xp2
        
        x1_next = x1_curr + 0.5 * (k1_x1 + k2_x1)
        x2_next = x2_curr + 0.5 * (k1_x2 + k2_x2)
        
        x1.append(x1_next)
        x2.append(x2_next)
        y.append(x2_curr)
        time.append(i * dT)
    
    return time, U[:-1], y
time_old, input_old, output_old = runge_kutta(generate_signal())

def update_params(val):
    
    """Update all parameters from sliders"""
    global L, C, R1, R2, amplitude, frequency, T, duty_cycle, input_old, output_old, time_old, K_old, param_old
    global harmonic_amplitudes, harmonic_frequencies
    global slider_L, slider_C, slider_R1, slider_R2, slider_amp, slider_freq, slider_time
    global slider_duty_cycle, slider_h1_amp, slider_h1_freq, slider_h2_amp, slider_h2_freq, slider_h3_amp, slider_h3_freq
    if(T == slider_time.val):
        time_old, input_old, output_old = runge_kutta(generate_signal())
    if(T==slider_time.val and frequency == slider_freq.val and amplitude == slider_amp.val and duty_cycle == slider_duty_cycle and harmonic_amplitudes[0] == slider_h1_amp and harmonic_amplitudes[1] == slider_h2_amp and harmonic_amplitudes[2] == slider_h3_amp and harmonic_frequencies[0] == slider_h1_freq and harmonic_frequencies[1] == slider_h2_freq and harmonic_frequencies[2] == slider_h3_freq):
        K_old = calc_coefficients()
        param_old = [L,C,R1,R2]
    # Update from sliders
    L = slider_L.val
    C = slider_C.val
    R1 = slider_R1.val
    R2 = slider_R2.val
    amplitude = slider_amp.val
    frequency = slider_freq.val
    T = slider_time.val
    duty_cycle = slider_duty_cycle.val
    
    harmonic_amplitudes[0] = slider_h1_amp.val
    harmonic_frequencies[0] = slider_h1_freq.val
    harmonic_amplitudes[1] = slider_h2_amp.val
    harmonic_frequencies[1] = slider_h2_freq.val  
    harmonic_amplitudes[2] = slider_h3_amp.val
    harmonic_frequencies[2] = slider_h3_freq.val
    
    update_plot()

def update_signal_type(label):
    
    """Update signal type from radio buttons"""
    global signal_type, input_old, output_old, time_old
    global slider_h1_amp, slider_h1_freq, slider_h2_amp, slider_h2_freq, slider_h3_amp, slider_h3_freq
    global slider_duty_cycle
    time_old, input_old, output_old = runge_kutta(generate_signal())
    
    signal_type = label
    
    # Pokaż/ukryj suwaki harmoniczne w zależności od wybranego typu sygnału
    if signal_type == 'harmonic':
        slider_h1_amp.ax.set_visible(True)
        slider_h1_freq.ax.set_visible(True)
        slider_h2_amp.ax.set_visible(True)
        slider_h2_freq.ax.set_visible(True)
        slider_h3_amp.ax.set_visible(True)
        slider_h3_freq.ax.set_visible(True)
    else:
        slider_h1_amp.ax.set_visible(False)
        slider_h1_freq.ax.set_visible(False)
        slider_h2_amp.ax.set_visible(False)
        slider_h2_freq.ax.set_visible(False)
        slider_h3_amp.ax.set_visible(False)
        slider_h3_freq.ax.set_visible(False)
    
    # Pokaż/ukryj suwak duty cycle dla sygnału prostokątnego
    if signal_type == 'square':
        slider_duty_cycle.ax.set_visible(True)
    else:
        slider_duty_cycle.ax.set_visible(False)
    
    plt.draw()
    update_plot()

def plot_frequency_response():
    """Plot magnitude response in the third subplot"""
    global ax3, K_old
    K = calc_coefficients()
    w = np.logspace(-1, 3, 500)
    
    magnitude_db = []
    magnitude_old = []
    for frequency in w:
        mag_old = 1 / (K_old[4] * np.sqrt(frequency**4 + K_old[5] * frequency**2 + K_old[6]))
        magnitude = 1 / (K[4] * np.sqrt(frequency**4 + K[5] * frequency**2 + K[6]))
        magnitude_old.append(20 * np.log10(mag_old))
        magnitude_db.append(20 * np.log10(magnitude))
    
    ax3.semilogx(w, magnitude_db, 'g-', linewidth=2, label='Magnitude')
    ax3.semilogx(w, magnitude_old, 'g--', alpha=0.5, label='Magnitude')
    ax3.set_xlabel('Frequency [rad/s]')
    ax3.set_ylabel('Magnitude [dB]')
    ax3.set_title('Magnitude Response')
    ax3.grid(True, alpha=0.3)

def plot_phase_response():
    """Plot phase response in the fourth subplot"""
    global ax4, param_old
    K = calc_coefficients()
    w = np.logspace(-1, 3, 500)
    
    phase_deg = []
    phase_old = []
    for frequency in w:
        # Phase: H(s) = 1 / (LCs^2 + (L/R1 + L/R2)s + 1)
        # Phase = -arctan(Im(H)/Re(H))
        s = 1j * frequency  # s = j*omega
        den_old = param_old[0] * param_old[1] * s**2 + (param_old[0]/param_old[2] + param_old[0]/param_old[3]) * s + 1
        denominator = L * C * s**2 + (L/R1 + L/R2) * s + 1
        phase_rad = -np.angle(denominator)
        phase_old.append(np.degrees(-np.angle(den_old)))
        phase_deg.append(np.degrees(phase_rad))
    
    ax4.semilogx(w, phase_deg, 'm-', linewidth=2, label='Phase')
    ax4.semilogx(w, phase_old, 'm--', alpha = 0.5, label='Phase')
    ax4.set_xlabel('Frequency [rad/s]')
    ax4.set_ylabel('Phase [degrees]')
    ax4.set_title('Phase Response')
    ax4.grid(True, alpha=0.3)

def update_plot():
    global ax1, ax2, ax3, ax4, controls_ax, L, C, R1, R2, amplitude, frequency, signal_type, input_old, output_old, time_old
    
    U = generate_signal()
    time, input_signal, output_signal = runge_kutta(U)
    
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    
    # input plot
    ax1.plot(time, input_signal, 'b-', linewidth=2, label='Input Signal')
    ax1.plot(time_old, input_old, 'b--', alpha = 0.5, label='Input Signal')
    ax1.set_ylabel('Amplitude [V]')
    ax1.set_title(f'Input Signal - {signal_type}')
    ax1.grid(True, alpha=0.3)
    
    # output plot
    ax2.plot(time, output_signal, 'r-', linewidth=2, label='Output Signal')
    ax2.plot(time_old, output_old, 'r--', alpha = 0.5, label='Output Signal')
    ax2.set_ylabel('Voltage [V]')
    ax2.set_title('Output Signal')
    ax2.grid(True, alpha=0.3)
    
    # freq plot
    plot_frequency_response()
    
    # phase plot
    plot_phase_response()
    
    # circuit info
    K = calc_coefficients()
    
    controls_ax.clear()
    controls_ax.set_xlim(0, 1)
    controls_ax.set_ylim(0, 1)
    controls_ax.axis('off')
    controls_ax.text(0, 1, "Interactive RLC Simulator", fontsize=15, 
                        verticalalignment='top')

    
    plt.draw()

fig = plt.figure(figsize=(18, 12))

gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3, width_ratios=[0.6, 0.5, 0.5]) # wspace=0.3, width_ratios=[0.6, 1, 1])


ax1 = fig.add_subplot(gs[0, 1])  # Input
ax2 = fig.add_subplot(gs[0, 2])  # Output
ax3 = fig.add_subplot(gs[1, 1])  # Magnitude
ax4 = fig.add_subplot(gs[1, 2])  # Phase


controls_ax = fig.add_subplot(gs[:, 0])
controls_ax.set_xlim(0, 1)
controls_ax.set_ylim(0, 1)
controls_ax.axis('off')

slider_props = dict(width=0.25, height=0.03)
slider_L = Slider(plt.axes([0.07, 0.75, slider_props['width'], slider_props['height']]), 'L [H]', 0.01, 10.0, valinit=L, valfmt='%.2f')
slider_C = Slider(plt.axes([0.07, 0.70, slider_props['width'], slider_props['height']]), 'C [F]', 1e-6, 1e-3, valinit=C, valfmt='%.2e')
slider_R1 = Slider(plt.axes([0.07, 0.65, slider_props['width'], slider_props['height']]), 'R1 [Ω]', 10, 10000, valinit=R1, valfmt='%.0f')
slider_R2 = Slider(plt.axes([0.07, 0.60, slider_props['width'], slider_props['height']]), 'R2 [Ω]', 10, 10000, valinit=R2, valfmt='%.0f')

slider_amp = Slider(plt.axes([0.07, 0.50, slider_props['width'], slider_props['height']]), 'Amplitude [V]', 0.1, 5.0, valinit=amplitude, valfmt='%.1f')
slider_freq = Slider(plt.axes([0.07, 0.45, slider_props['width'], slider_props['height']]), 'Freq [Hz]', 1, 100, valinit=frequency, valfmt='%.0f')
slider_time = Slider(plt.axes([0.07, 0.40, slider_props['width'], slider_props['height']]), 'Time [s]', 0.05, 1.0, valinit=T, valfmt='%.2f')

slider_duty_cycle = Slider(plt.axes([0.07, 0.35, slider_props['width'], slider_props['height']]), 'Duty Cycle [%]', 0, 100, valinit=duty_cycle, valfmt='%.0f')
slider_duty_cycle.ax.set_visible(False)


slider_h1_amp = Slider(plt.axes([0.07, 0.32, slider_props['width'], slider_props['height']]), 'H1 Amp', 0.0, 2.0, valinit=harmonic_amplitudes[0], valfmt='%.1f')

slider_h1_freq = Slider(plt.axes([0.07, 0.29, slider_props['width'], slider_props['height']]), 'H1 Freq [Hz]', 1, 100, valinit=harmonic_frequencies[0], valfmt='%.0f')

slider_h2_amp = Slider(plt.axes([0.07, 0.26, slider_props['width'], slider_props['height']]), 'H2 Amp', 0.0, 2.0, valinit=harmonic_amplitudes[1], valfmt='%.1f')

slider_h2_freq = Slider(plt.axes([0.07, 0.23, slider_props['width'], slider_props['height']]), 'H2 Freq [Hz]', 1, 100, valinit=harmonic_frequencies[1], valfmt='%.0f')

slider_h3_amp = Slider(plt.axes([0.07, 0.20, slider_props['width'], slider_props['height']]), 'H3 Amp', 0.0, 2.0, valinit=harmonic_amplitudes[2], valfmt='%.1f')

slider_h3_freq = Slider(plt.axes([0.07, 0.17, slider_props['width'], slider_props['height']]), 'H3 Freq [Hz]', 1, 100, valinit=harmonic_frequencies[2], valfmt='%.0f')

# Ukryj suwaki harmoniczne domyślnie (ponieważ domyślnie wybrana jest opcja 'step')
slider_h1_amp.ax.set_visible(False)
slider_h1_freq.ax.set_visible(False)
slider_h2_amp.ax.set_visible(False)
slider_h2_freq.ax.set_visible(False)
slider_h3_amp.ax.set_visible(False)
slider_h3_freq.ax.set_visible(False)

# Signal type radio buttons
radio_signal = RadioButtons(plt.axes([0.07, 0.02, 0.15, 0.12]), ('pulse', 'step', 'sine', 'square', 'saw', 'harmonic'),active=1)

# Connect slider events
slider_L.on_changed(update_params)
slider_C.on_changed(update_params)
slider_R1.on_changed(update_params)
slider_R2.on_changed(update_params)
slider_amp.on_changed(update_params)
slider_freq.on_changed(update_params)
slider_time.on_changed(update_params)
slider_duty_cycle.on_changed(update_params)

# Connect harmonic sliders
slider_h1_amp.on_changed(update_params)
slider_h1_freq.on_changed(update_params)
slider_h2_amp.on_changed(update_params)
slider_h2_freq.on_changed(update_params)
slider_h3_amp.on_changed(update_params)
slider_h3_freq.on_changed(update_params)
radio_signal.on_clicked(update_signal_type)

update_plot()

plt.show()