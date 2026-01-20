## Zadání
![Zadání ulohy](/images/uloha1.png "Zadání úlohy")

## Kód
```r 
# ==============================================================================
# 1. DEFINICE FUNKCÍ (Numerické metody)
# ==============================================================================

# Gaussova eliminace s pivotací
GaussEliminationPivoting <- function(A, b) {
  N <- length(b)
  Ab <- cbind(A, b)  # rozšířená matice [A | b]
  
  for (p in 1:(N-1)) {
    # Pivotace
    imax <- which.max(abs(Ab[p:N, p])) + p - 1
    if (abs(Ab[imax, p]) < .Machine$double.eps) stop("Singulární matice.")
    
    if (imax != p) {
      Ab[c(p, imax), ] <- Ab[c(imax, p), ] # Prohození řádků
    }
    
    # Eliminace
    for (r in (p+1):N) {
      m <- Ab[r, p] / Ab[p, p]
      Ab[r, p:(N+1)] <- Ab[r, p:(N+1)] - m * Ab[p, p:(N+1)]
    }
  }
  
  # Zpětný chod
  x <- numeric(N)
  x[N] <- Ab[N, N+1] / Ab[N, N]
  for (r in (N-1):1) {
    x[r] <- (Ab[r, N+1] - sum(Ab[r, (r+1):N] * x[(r+1):N])) / Ab[r, r]
  }
  return(x)
}

# Lagrangeova interpolace
# Funkce Lagrange – vrací hodnotu interpolačního polynomu v bodě t
# Parametry:
#   t ... bod, kde chceme polynom vyhodnotit
#   x ... vektor uzlů (zadané hodnoty x_i)
#   y ... vektor odpovídajících hodnot funkce y_i

Lagrange <- function(t, x, y){
  n <- length(x)
  soucet <- 0
  for(i in 1:n){
    soucin <- 1
    for(j in 1:n){
      if(j != i) {
        soucin <- soucin * (t - x[j]) / (x[i] - x[j])
      }
    }
    soucet <- soucet + y[i] * soucin
  }
  return(soucet)
}

# Simpsonovo pravidlo
simpson <- function(f, a, b, n = 1) {
  h <- (b - a) / n
  # Simpson vyžaduje sudé n, pokud ne, R funkce integrate je robustnější, 
  # ale pro školní účely zde necháváme tvou implementaci.
  # Vzorec: (h/3) * (f(a) + f(b) + 4*sum(liché) + 2*sum(sudé))
  # Tvá implementace je mírně nestandardní, ale matematicky ekvivalentní pro složené pravidlo:
  suma <- f(a) + f(b) + 4 * sum(f(a + h * (1:n) - h / 2))
  if (n > 1) {
    suma <- suma + 2 * sum(f(a + h * (1:(n - 1))))
  }
  return(h * suma / 6)
}

# ==============================================================================
# 2. ŘEŠENÍ ÚLOHY
# ==============================================================================

# --- Zadání dat ---
n <- 10
A <- matrix(0, n, n)
b <- numeric(n) # Toto je vektor y ze zadání rovnic (3)

for (i in 1:n) {
  b[i] <- sin(i)
  for (j in 1:n) {
    A[i,j] <- cos((i-1)*j) - j
  }
}

# --- Úkol 1: Řešení soustavy Ax = y ---
# V zadání je vektor neznámých označen jako 'x'
vec_x <- GaussEliminationPivoting(A, b)

# Pro interpolaci máme body [x, y]:
# x-ové souřadnice jsou řešením soustavy (vec_x)
# y-ové souřadnice jsou pravou stranou soustavy (b)
x_nodes <- vec_x
y_nodes <- b

cat("Vypočítané řešení x (uzly interpolace):\n")
print(x_nodes)

# --- Úkol 2a: Průsečíky a Graf ---

# Definujeme funkci polynomu p(t) s fixovanými uzly pro další výpočty
poly_p <- function(t) {
  Lagrange(t, x_nodes, y_nodes)
}

# 1. Průsečík s osou Y (hodnota v bodě 0)
prusecik_Y <- poly_p(0)
cat(sprintf("\nPrůsečík s osou Y [0, p(0)]: [0, %.6f]\n", prusecik_Y))


# --- Úkol 2b: Integrál ---
# Meze integrálu podle zadání (min a max prvek vektoru x)
x_min <- min(x_nodes)
x_max <- max(x_nodes)

# Výpočet integrálu (Wrapper 'poly_p' musí být vektorizovaný pro Simpsonovo pravidlo uvnitř)
integrand_vec <- Vectorize(poly_p)

# Použijeme dostatečný počet kroků (n=200)
integral_result <- simpson(integrand_vec, x_min, x_max, n = 200)

cat(sprintf("\n--- Integrál ---\n"))
cat(sprintf("Meze: [%.6f, %.6f]\n", x_min, x_max))
cat(sprintf("Hodnota integrálu: %.6f\n", integral_result))


# ==============================================================================
# 3. VYKRESLENÍ GRAFUhttp://127.0.0.1:39215/graphics/plot_zoom_png?width=1201&height=583
# ==============================================================================

# Generování bodů pro hladkou křivku
tt <- seq(min(x_nodes) - 1, max(x_nodes) + 1, length.out = 300)
pp <- sapply(tt, poly_p)

# Graf
plot(x_nodes, y_nodes, col = "blue", pch = 16, cex = 1.2,
     xlim = c(-5,5), 
     ylim = c(-15,20), 
     xlab = "x", ylab = "p(x)",
     main = "Lagrangeův interpolační polynom")

# Křivka polynomu
lines(tt, pp, col = "red", lwd = 2)

# Mřížka a osy
grid()
abline(h = 0, v = 0, col = "gray50", lty = 2)

# Vyznačení průsečíku s osou Y (zeleně)
points(0, prusecik_Y, col = "green", pch = 19, cex = 1.5)
text(0, prusecik_Y, labels = sprintf("  Y-intersect\n  %.2f", prusecik_Y), adj = c(0, 1), col="green4")

# Vyznačení plochy integrálu (stínování)
polygon_x <- c(tt[tt >= x_min & tt <= x_max], x_max, x_min)
polygon_y <- c(pp[tt >= x_min & tt <= x_max], 0, 0)
polygon(polygon_x, polygon_y, col = rgb(1, 0, 0, 0.1), border = NA)
```
## Výstup

![Výstup v R](/images/vystup1.png "Zadání úlohy")