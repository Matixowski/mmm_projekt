# Wnioski z pracy projektowej

## Działanie algorytmu Rungego-Kutty

Algorytm Rungego-Kutty 2. rzędu pozwala na rozwiązanie równań różniczkowych potrzebnych do realizacji układu. Został wybrany ze względu na swoją dokładność w porównaniu do metod Eulera, Heuena i innych [[porównanie z Wikipedii]](https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Runge-kutta.svg/1920px-Runge-kutta.svg.png) oraz prostotę w zastosowaniu.

### Zasada działania algorytmu RK2

Algorytm Rungego-Kutty 2. rzędu rozwiązuje równania stanu obwodu RLC:
- **x1**: prąd przez cewkę (I_L)  
- **x2**: napięcie na kondensatorze (U_C)

**Procedura obliczeń w każdym kroku:**

1. **Obliczenie współczynników k1**:
   - k1_x1 = K[0] * x2_curr + K[3] * u_curr
   - k1_x2 = K[1] * x1_curr + K[2] * x2_curr

2. **Przewidywane wartości**:
   - xp1 = x1_curr + k1_x1
   - xp2 = x2_curr + k1_x2

3. **Obliczenie współczynników k2**:
   - k2_x1 = K[0] * xp2 + K[3] * u_next
   - k2_x2 = K[1] * xp1 + K[2] * xp2

4. **Końcowa aktualizacja stanu**:
   - x1_next = x1_curr + 0.5 * (k1_x1 + k2_x1)
   - x2_next = x2_curr + 0.5 * (k1_x2 + k2_x2)

### Współczynniki macierzy stanu

Algorytm wykorzystuje prekalkulowane współczynniki K dla reprezentacji stanu:
- **K[0] = -dT/L**: wpływ napięcia kondensatora na prąd cewki
- **K[1] = dT/C**: wpływ prądu cewki na napięcie kondensatora  
- **K[2] = -dT/(R1*C) - dT/(R2*C)**: rozładowanie kondensatora przez rezystory
- **K[3] = dT/L**: wpływ sygnału wejściowego na prąd cewki




## Wpływ zmiany parametrów sygnału wejściowego na:

### Napięcie wyjściowe, wzmocnienie i przesunięcie fazowe

**Parametry sygnału:**
- **Amplituda:** Wpływa bezpośrednio na amplitudę wyjściową (z wyjątkiem sygnału harmonicznego, gdzie znaczenie mają amplitudy poszczególnych składowych). Nie ma wpływu na charakterystyki wzmocnienia ani przesunięcia fazowego.
- **Częstotliwość i czas:** Dla sygnałów zmiennych pozwalają zaobserwować ustalenie się sygnału przed kolejną zmianą (np. impulsem sygnału prostokątnego). Zwiększenie czasu symulacji umożliwia obserwację sygnału ustalonego. Nie wpływają na charakterystyki wzmocnienia ani przesunięcia fazowego.

**Parametry obwodu:**
- **Rezystancja (R1, R2):** Zwiększenie oporu powoduje nieliniowy wzrost napięcia wyjściowego, podniesienie punktu częstotliwości rezonansowej na charakterystyce wzmocnienia oraz zwiększenie "stromizny" przesunięcia fazowego.
- **Indukcyjność (L):** Zwiększenie indukcyjności skraca czas stabilizacji sygnału wyjściowego, zwiększa dobroć układu oraz wypłaszcza "stromiznę" fazy, czyniąc ją bardziej liniową.
- **Pojemność (C):** Zwiększenie pojemności wydłuża czas stabilizacji sygnału wyjściowego, obniża punkt częstotliwości rezonansowej oraz zmniejsza moment odwrócenia fazy.

Dla wysokich częstotliwości przesunięcie fazowe dąży do -180°.




## Analiza typów sygnałów wejściowych

### Sygnał impulsowy (pulse)
- **Charakterystyka**: Pojedynczy impuls o amplitudzie A na początku symulacji
- **Odpowiedź**: Oscylacje wygaszane do zera

### Sygnał skokowy (step)
- **Charakterystyka**: Stały sygnał o amplitudzie A przez cały czas symulacji
- **Odpowiedź**: Oscylacje wygaszane do wartości ustalonej

### Sygnał sinusoidalny (sine)
- **Charakterystyka**: u(t) = A * sin(2πft)
- **Odpowiedź**: Sygnał sinusoidalny z wygaszanymi oscylacjami, po ustabilizowaniu przyjmuje postać sinusoidy o tej samej częstotliwości co wejściowa

### Sygnał prostokątny (square)
- **Charakterystyka**: Przełączanie między amplitudą A a 0V z częstotliwością f, przy czym stosunek czasu wysokiego do okresu określa duty_cycle
- **Parametr dodatkowy**: duty_cycle [%]
- **Odpowiedź**: Oscylacje przy każdej zmianie poziomu sygnału, wygaszane do momentu kolejnej zmiany

### Sygnał piłokształtny (saw)
- **Charakterystyka**: u(t) = A * (2*(t%T)/T - 1)
- **Odpowiedź**: Liniowy wzrost i nagły spadek napięcia wyjściowego z nałożonymi oscylacjami, które są wygaszane w trakcie liniowego narastania sygnału a następnie znowu wracają po zmianie

### Sygnał harmoniczny (harmonic)
- **Charakterystyka**: u(t) = A * suma[H_i * sin(2* pi * f_i * t + fi_i)] dla i=1,2,3
- **Parametry dodatkowe**: amplitude [V] oraz freq [Hz] dla każdego z trzech sygnałów składowych
- **Odpowiedź**: Złożony sygnał wyjściowy będący kombinacją odpowiedzi na poszczególne składowe harmoniczne, z różnymi wzmocnieniami i przesunięciami fazowymi zależnymi od częstotliwości każdej składowej
