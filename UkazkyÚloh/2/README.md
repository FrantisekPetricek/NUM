## Zadání

![Zadání ulohy](/images/uloha2.png "Zadání úlohy")

## Kód

```r
midpointrule <- function(f, a, b, n = 100) { 
  h <- (b - a) / n
  return(h * sum(f(a + h * (1:n) - 0.5 * h)))
}

steffensen <- function(f, x0, tol = 1e-8, max_iter = 100, verbose = TRUE) {
  if (verbose) {
    cat("Iterace |        x         |        f(x)        \n")
    cat("------------------------------------------------\n")
  }
  for (i in 1:max_iter) {
    fx <- f(x0)
    if(!is.finite(fx)) stop("Funkce vrátila nekonečno nebo NaN.")
    
    h  <- fx                       # odhad kroku (h = f(x))
    
    # Ochrana proti dělení nulou, pokud jsme už u kořene
    if (abs(h) < 1e-14) {
      return(x0)
    }
    
    gx <- (f(x0 + h) - fx) / h     # g(x) ≈ kvaziderivace
    
    x1 <- x0 - fx / gx             # Steffensenův krok
    
    if (verbose) {
      cat(sprintf("%5d   | %14.10f | %14.6e\n", i, x1, fx))
    }
    
    if (abs(x1 - x0) < tol) {
      if (verbose) {
        cat("\n Konvergence dosažena po", i, "iteracích.\n")
        cat(sprintf("Kořen: x = %.12f\n", x1))
      }
      return(x1)
    }
    x0 <- x1
  }
  return(x1)
}

# ---------------------------------------------------------
# 2. DEFINICE PROBLÉMU PRO NUMERICKÉ ŘEŠENÍ
# ---------------------------------------------------------

# Cílová hodnota ze zadání
TARGET_VAL <- 1.31848

# A) Definice vnitřního integrálu (podle x)
# Vypočítá hodnotu I(alpha) = integral_0^2pi (e^(-alpha*x) * sin(x)) dx
calc_inner_integral <- function(alpha) {
  
  # Integrand uvnitř: funkce proměnné x, alpha je parametr
  integrand_x <- function(x) {
    exp(-alpha * x) * sin(x)
  }
  
  # Použijeme tvůj midpoint rule. 
  # n=200 pro dostatečnou přesnost na vlně sinusu
  midpointrule(integrand_x, a = 0, b = 2*pi, n = 200)
}

# DŮLEŽITÉ: Tvůj 'midpointrule' posílá do funkce 'f' vektor hodnot (x body).
# Funkce 'calc_inner_integral' ale umí zpracovat jen jedno alpha najednou.
# Abychom ji mohli použít jako 'f' ve vnějším integrálu, musíme ji "vektorizovat".
calc_inner_integral_vec <- Vectorize(calc_inner_integral)


# B) Definice vnějšího integrálu (podle alpha)
# Vypočítá integral_0^p ( I(alpha) ) d_alpha
calc_total_integral <- function(p) {
  # Pokud je p blízko 0, integrál je 0
  if(abs(p) < 1e-9) return(0)
  
  # Použijeme tvůj midpoint rule na vektorizovanou vnitřní funkci
  # n=200 opět pro přesnost
  midpointrule(calc_inner_integral_vec, a = 0, b = p, n = 200)
}

# C) Cílová funkce pro Steffensenovu metodu
# Hledáme p, kde f_to_solve(p) = 0
f_to_solve <- function(p) {
  # Spočítáme dvojný integrál numericky a odečteme cílovou hodnotu
  val <- calc_total_integral(p)
  return(val - TARGET_VAL)
}

# ---------------------------------------------------------
# 3. SPUŠTĚNÍ VÝPOČTU
# ---------------------------------------------------------

cat("\n--- Hledání p pro dvojný integrál ---\n")
# Zvolíme počáteční odhad. Zkusíme například p = 2.
# (Víme, že řešení je kolem 10, ale nechme metodu pracovat)
vysledek_p <- steffensen(f_to_solve, x0 = 2.0, tol = 1e-6)

cat(sprintf("\nVýsledné p: %.5f\n", vysledek_p))

# Rychlá kontrola výsledku zpětným dosazením
cat(sprintf("Hodnota integrálu pro nalezené p: %.6f (Cíl: %.6f)\n", 
            calc_total_integral(vysledek_p), TARGET_VAL))

```

## Výstup


```text
--- Hledání p pro dvojný integrál ---
Iterace |        x         |        f(x)        
------------------------------------------------
    1   |   3.5579002868 |  -3.638734e-01
    2   |   5.8256225878 |  -1.741011e-01
    3   |   8.2286869390 |  -6.987203e-02
    4   |   9.6188510132 |  -2.048759e-02
    5   |   9.8882685205 |  -2.926403e-03
    6   |   9.8958536875 |  -7.807708e-05
    7   |   9.8958593784 |  -5.849142e-08
    8   |   9.8958593784 |  -2.597922e-14

    Konvergence dosažena po 8 iteracích.
Kořen: x = 9.895859378379

Výsledné p: 9.89586

Hodnota integrálu pro nalezené p: 1.318480 (Cíl: 1.318480)
```