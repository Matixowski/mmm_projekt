# Projekt MMM
## Autorzy
* Mateusz Gniado (197888)
* Piotr Szulc (197639)
## Treść zadania (Projekt nr 16)

Dany jest układ RLC

 ```
--- L ------------
 |       |   |    |
 |       |   |    |
 u       R   C    R <- y(t) = U_R
 |       |   |    |
 ------------------
 ```

gdzie R, L, C to parametry modelu a u to sygnał wejściowy. Należy zaimplementować symulator tego układu
umożliwiając uzyskanie odpowiedzi czasowych układu na pobudzenie przynajmniej trzema
rodzajami synagłów wejściowych (prostokątny o skończonym czasie trwania, trójkątny,
harmoniczny). Symulator powinien umożliwiać zmianę wszystkich parametrów modelu oraz
sygnałów wejściowych. Należy wykreślić charakterystyki częstotliwościowe Bodego
(amplitudową i fazową) oraz odpowiedź układu.

## Jak używać

Po uruchomieniu wyświetli się GUI MatPlota. Przesuwając suwakami można ustawiać parametry programu. Również można przyciskami zmienić typ sygnału na wejściu. 
Przy włączeniu sygnału prostokątnego wyświetli się dodatkowy suwak "duty cycle" a przy włączeniu harmonicznego będzie możliwość zmiany amplitudy i częstotliwości dla trzech składowych.

## Instalacja

1. Pobierz projekt
2. Zainstaluj matplotlib za pomocą `pip install matplotlib`
3. Uruchom za pomocą `python rlc_simulator.py`

## Wnioski z pracy projektowej

Wnioski są dostępne w pliku [wnioski.md](wnioski.md).

## Architektura i implementacja programu

### Funkcje
1. **generate_signal()**: Generacja różnych typów sygnałów wejściowych
2. **calc_coefficients()**: Obliczanie współczynników macierzy stanu
3. **runge_kutta()**: Numeryczne rozwiązanie równań różniczkowych (wymagane w realizacji projektu)
4. **plot_frequency_response()** / **plot_phase_response()**: Analiza częstotliwościowa
5. **update_params()** / **update_signal_type()**: Obsługa GUI

### Interface użytkownika
**Organizacja kontrolek:**
- **Parametry obwodu**: L, C, R1, R2
- **Parametry sygnału**: Amplituda, częstotliwość, czas
- **Parametry dodatkowe**: duty_cycle (square), harmoniczne [po amplitudzie i częstotliwości dla każdej z trzech składowych] (harmonic)
- **Typ sygnału**: Wybór sześciu rodzajów sygnałów
- **Dynamiczna widoczność**: Automatyczne pokazywanie/ukrywanie odpowiednich suwaków

### Wizualizacja wyników
**Układ 4 wykresów:**
1. **Sygnał wejściowy**: (gs[0,1])
2. **Sygnał wyjściowy**: (gs[0,2])
3. **Charakterystyka amplitudowa**: (gs[1,1])
4. **Charakterystyka fazowa**: (gs[1,2])