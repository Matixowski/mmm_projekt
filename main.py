import matplotlib.pyplot as plt
import numpy as np

def calcCoeff(x, dT):
    a = -1*dT/x[0]
    b = 1*dT/x[1]
    c = -1*dT/(x[2]*x[1]) - 1*dT/(x[3]*x[1])
    d = 1*dT/x[0]
    g1 = x[0]*x[1]
    g2 = (x[0]*(x[2]+x[3])**2 - 2*(x[2]**2)*(x[3]**2)*x[1])/((x[2]**2)*(x[3]**2)*(x[1]**2)*x[0])
    g3 = 1/((x[0]**2)*(x[1]**2))
    p1 = 1/(x[2]*x[1]) + 1/(x[3]*x[1])
    p2 = 1/(x[0]*x[1])
    condition = np.sqrt(1/(x[0]*x[1]))
    return [a, b, c, d, g1, g2, g3, p1, p2, condition]

def prepInput(T, dT):
    U = []
    i = 0
    while i < T/dT:
        U.append(1)
        i += 1
    return U

# A = [0 a]
#     [b c]

# B = [d]
#     [0]

#|G(jw)| = 1/(g1*sqrt(w^4 + g2*w^2 + g3))

#arg(G(jw)) = arctg((p1*w)/(w^2 - p2))

#    a  b  c  d g1 g2 g3 p1 p2  
K = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

T = 0.2
dT = 0.001

#    L     C      R1    R2
P = [1, 0.00001, 1000, 1000]
K = calcCoeff(P, dT)
U = prepInput(T,dT)

x1 = [0]
xp1 = []
x2 = [0]
xp2 = []
y = []
time = []

i = 0

gain = []
phase = []
puls = []
ph = 0
w = 0

while i<(T/dT)-1:
    xp1.append(x1[i] + K[0]*x2[i] + K[3]*U[i])
    xp2.append(x2[i] + K[1]*x1[i] + K[2]*x2[i])
    x1.append(x1[i] + 0.5*(K[0]*x2[i] + K[3]*U[i] + K[0]*xp2[i] + K[3]*U[i+1]))
    x2.append(x2[i] + 0.5*(K[1]*x1[i] + K[2]*x2[i] + K[1]*xp1[i] + K[2]*xp2[i]))
    y.append(x2[i])
    time.append(i*dT)
    i += 1
while ph>-179:
    gain.append(20*np.log10(1/(K[4]*np.sqrt(w**4+K[5]*w**2+K[6]))))
    ph = 180*np.arctan(K[7]*w/(w**2-K[8]))/np.pi
    if w>K[9]:
        ph -= 180
    phase.append(ph)
    puls.append(w)
    w += 0.1

f1 = plt.figure(1)
plt.plot(time, y)
plt.xlabel("Czas [s]")
plt.ylabel("Napięcie [V]")
f2 = plt.figure(2)
plt.plot(puls, phase)
plt.xscale('log')
plt.xlabel("Pulsacja [rad/s]")
plt.ylabel("Przesunięcie fazowe [°]")
f3 = plt.figure(3)
plt.plot(puls, gain)
plt.xscale('log')
plt.xlabel("Pulsacja [rad/s]")
plt.ylabel("Wzmocnienie [dB]")
plt.show()