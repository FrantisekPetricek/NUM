# Detailní přehled numerických metod

## 1. Metody hledání kořene (Root finding)
Cílem je nalézt číslo $x$, pro které platí $f(x) = 0$. Geometricky hledáme průsečík křivky s vodorovnou osou x.

### Metoda půlení intervalu (Bisection method)
**Princip:**
Tato metoda je založena na Bolzanově větě. Pokud je funkce spojitá na intervalu $[a, b]$ a platí $f(a) \cdot f(b) < 0$ (funkce má na koncích opačná znaménka), musí uvnitř existovat kořen.
1. Najdeme střed intervalu $c = (a+b)/2$.
2. Zjistíme znaménko $f(c)$.
3. Pokud má $f(a)$ a $f(c)$ opačná znaménka, kořen je v levé polovině. Pokud ne, je v pravé.
4. Proces opakujeme, dokud není interval dostatečně malý.
**Výhody/Nevýhody:** Vždy konverguje (je robustní), ale je velmi pomalá.


### Metoda tečen (Newtonova metoda / Newton-Raphson)
**Princip:**
Místo půlení intervalu využíváme derivaci funkce (směrnici).
1. Zvolíme počáteční odhad $x_0$.
2. V tomto bodě sestrojíme **tečnu** ke grafu funkce.
3. Průsečík tečny s osou x je náš nový, lepší odhad $x_1$.
4. Vzorec: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$.
**Výhody/Nevýhody:** Konverguje extrémně rychle (kvadraticky) v blízkosti kořene. Vyžaduje však výpočet derivace a může selhat, pokud je derivace nulová nebo startovní bod špatný.


### Steffensenova metoda
**Princip:**
Jedná se o metodu pro hledání pevného bodu ($x = g(x)$), která je vylepšena o tzv. **Aitkenův $\Delta^2$ proces**.
Metoda se snaží "předpovědět", kam by iterace doputovala v budoucnu, kdybychom pokračovali v prostém dosazování, a rovnou skočí na tento odhad.
**Výhody/Nevýhody:** Dosahuje stejné rychlosti jako Newtonova metoda (kvadratická konvergence), ale **nevyžaduje výpočet derivace**, což je obrovská výhoda u složitých funkcí.


### Halleyova metoda
**Princip:**
Je to rozšíření Newtonovy metody. Zatímco Newton používá lineární aproximaci (tečnu), Halleyova metoda využívá i **druhou derivaci** a aproximuje funkci hyperbolou. Bere tedy v úvahu i zakřivení (konvexitu/konkavitu) funkce.
**Výhody/Nevýhody:** Má kubickou konvergenci (ještě rychlejší než Newton), ale vyžaduje výpočet první i druhé derivace, což je výpočetně náročné.


### Regula falsi (Metoda sečen)
**Princip:**
Vylepšení metody půlení intervalu. Máme dva body $a, b$ s opačnými znaménky funkčních hodnot.
Místo abychom interval rozpůlili v polovině, spojíme body $[a, f(a)]$ a $[b, f(b)]$ přímkou (**sečnou**). Nový bod leží tam, kde sečna protne osu x.
**Výhody/Nevýhody:** Zpravidla rychlejší než půlení intervalu, vždy konverguje, ale může být pomalá, pokud je funkce "prohnutá" pouze na jednu stranu.


### Fixed point iteration (Metoda prosté iterace)
**Princip:**
Rovnici $f(x)=0$ přepíšeme na tvar $x = g(x)$. Hledáme průsečík funkce $g(x)$ a přímky $y=x$.
Začneme s $x_0$, vypočteme $x_1 = g(x_0)$, pak $x_2 = g(x_1)$ atd.
**Výhody/Nevýhody:** Velmi jednoduchá na implementaci. Konverguje pouze tehdy, pokud je derivace $|g'(x)| < 1$ v okolí kořene (funkce nesmí růst příliš strmě). Graficky se znázorňuje jako "pavučinový graf".


### Newtonova-Hornerova metoda
**Princip:**
Specifická varianta Newtonovy metody určená výhradně pro hledání kořenů **polynomů**.
K vyčíslení hodnoty polynomu $P(x)$ a jeho derivace $P'(x)$ v daném bodě se používá **Hornerovo schéma**. To je algoritmus, který minimalizuje počet násobení a je numericky stabilnější než prosté dosazování do mocnin.


---

## 2. Řešení soustav lineárních rovnic
Hledáme vektor $\vec{x}$, aby platilo $A\vec{x} = \vec{b}$.

### Gaussova eliminační metoda (GEM)
**Princip:**
1. **Přímý chod:** Pomocí ekvivalentních úprav (prohazování řádků, přičítání násobků) upravíme matici na **horní trojúhelníkový tvar** (pod hlavní diagonálou jsou samé nuly).
2. **Zpětný chod:** Od posledního řádku (kde je rovnice s jednou neznámou) vypočítáme $x_n$, dosadíme do předposledního, vypočítáme $x_{n-1}$ atd.


### Metoda LU rozkladu
**Princip:**
Matici $A$ rozložíme na součin dvou matic: $A = L \cdot U$.
* $L$ (Lower) je dolní trojúhelníková matice (nuly nad diagonálou).
* $U$ (Upper) je horní trojúhelníková matice (nuly pod diagonálou).
Soustava $Ax=b$ se pak řeší nadvakrát: $Ly=b$ (hledáme y) a následně $Ux=y$ (hledáme x).
**Výhoda:** Pokud řešíme soustavu opakovaně pro stejnou matici $A$ ale různé vektory $b$, uděláme rozklad jen jednou a ušetříme čas.


### Jacobiho a Gaussova-Seidelova metoda
**Princip:**
Jedná se o **iterační metody** (postupně zpřesňují odhad, nepočítají přesně jako GEM). Z $i$-té rovnice vyjádříme $x_i$.
* **Jacobi:** Pro výpočet nové generace $x^{(k+1)}$ používáme na pravé straně pouze hodnoty ze staré generace $x^{(k)}$.
* **Gauss-Seidel:** Jakmile v cyklu vypočítáme nové $x_1$, okamžitě ho použijeme pro výpočet $x_2$ v téže iteraci.
**Výhoda:** Gauss-Seidel konverguje obvykle rychleji. Metody jsou skvělé pro velké řídké matice.


---

## 3. Interpolace
Máme diskrétní body (data) a hledáme funkci, která jimi **přesně prochází**.

### Lineární interpolace
**Princip:**
Nejjednodušší metoda. Sousední body spojíme úsečkami (přímkami). Vznikne lomená čára.
Používá se, pokud máme hustou síť bodů a nepotřebujeme hladkou křivku.


### Interpolace polynomem v Lagrangeově tvaru
**Princip:**
Pro $n+1$ bodů existuje právě jeden polynom stupně $n$, který jimi prochází.
Lagrangeův tvar tento polynom konstruuje jako vážený součet bázových polynomů $L_i(x)$.
Vlastnost $L_i(x)$ je, že v bodě $x_i$ má hodnotu 1 a ve všech ostatních uzlech 0.
**Nevýhoda:** Přidání jednoho nového bodu vyžaduje přepočítání celého polynomu.


### Interpolace polynomem v Newtonově tvaru
**Princip:**
Výsledný polynom je matematicky totožný s Lagrangeovým (stejná křivka), ale je zapsán v jiném tvaru:
$P(x) = a_0 + a_1(x-x_0) + a_2(x-x_0)(x-x_1) + \dots$
Koeficienty $a_i$ se počítají pomocí tabulky **poměrných diferencí**.
**Výhoda:** Pokud přidáme nový bod, stačí dopočítat jeden člen na konec, není třeba začínat od nuly.


---

## 4. Aproximace
Máme "mrak" bodů (např. z měření s chybou) a nechceme je spojit (graf by byl zubatý), ale proložit jimi hladkou křivku, která vystihuje **trend**.

### Metoda nejmenších čtverců (Least Squares)
**Princip:**
Hledáme parametry funkce (např. přímky $y=ax+b$), aby součet druhých mocnin (čtverců) svislých vzdáleností (reziduí) mezi naměřenými body a křivkou byl minimální.
Geometricky: Minimalizujeme plochu čtverečků mezi body a přímkou.


---

## 5. Derivace

### Taylorův rozvoj
**Princip:**
Každou "hezkou" (hladkou) funkci lze v okolí bodu nahradit nekonečným polynomem (mocninnou řadou), jehož koeficienty závisí na derivacích funkce.
$f(x) \approx f(a) + f'(a)(x-a) + \frac{f''(a)}{2!}(x-a)^2 + \dots$
V numerice se využívá k odvození vzorců pro numerickou derivaci (např. centrální diference) a k analýze chyb metod.


---

## 6. Integrace (Numerická kvadratura)
Počítáme určitý integrál $\int_a^b f(x) dx$, což je **obsah plochy pod grafem**.

### Midpoint rule (Obdélníková metoda)
**Princip:**
Interval rozdělíme na podintervaly. V každém z nich aproximujeme funkci konstantou – hodnotou v **půli** intervalu.
Plocha pod grafem je součtem ploch těchto obdélníků.


### Lichoběžníkové pravidlo (Trapezoidal rule)
**Princip:**
Interval rozdělíme na podintervaly. V každém z nich nahradíme funkci **úsečkou** spojující funkční hodnoty na krajích.
Plocha se počítá jako součet ploch lichoběžníků. Je přesnější než obdélníková metoda pro funkce, které nejsou konstantní.


### Simpsonova metoda
**Princip:**
Místo úseček (lineární) používáme proloženou **parabolu** (kvadratická funkce) vždy přes tři sousední body.
Vzorec využívá váhy 1-4-1. Je výrazně přesnější než lichoběžníkové pravidlo pro hladké funkce.


### Integrace pomocí Gaussovy kvadratury
**Princip:**
Zatímco předchozí metody (Newton-Cotesovy) mají body rozmístěné rovnoměrně, Gaussova kvadratura si vybírá "nejlepší" body (uzly) a váhy.
Uzly jsou kořeny Legendreových polynomů. Díky tomu dokáže jen s pár body přesně integrovat polynomy vysokých stupňů.


### Integrace pomocí Rombergovy kvadratury
**Princip:**
Chytrá metoda, která kombinuje **lichoběžníkové pravidlo** s **Richardsonovou extrapolací**.
Spočítá integrál hrubě (velký krok), pak jemněji (poloviční krok) a matematicky odhadne, k jaké hodnotě by se došlo, kdyby byl krok nulový. Tím eliminuje chybu a velmi rychle zpřesňuje výsledek.


---

## 7. Diferenciální rovnice (ODE)
Řešíme počáteční úlohu: známe stav v čase $t_0$ a rovnici pro rychlost změny $y' = f(t,y)$. Chceme predikovat budoucnost.

### Eulerova metoda
**Princip:**
Základní metoda. V aktuálním bodě spočítáme derivaci (sklon tečny) a uděláme malý krok $h$ po této přímce.
$y_{n+1} = y_n + h \cdot f(t_n, y_n)$.
Je málo přesná, protože ignoruje změnu sklonu během kroku.


### RK4 (Runge-Kutta 4. řádu)
**Princip:**
Standard v simulacích. V rámci jednoho časového kroku se "podívá" na sklon ve 4 místech:
1. Na začátku intervalu ($k_1$).
2. Dvakrát uprostřed intervalu (odhadnuto pomocí $k_1$ a pak $k_2$).
3. Na konci intervalu ($k_4$).
Výsledný krok je váženým průměrem těchto sklonů. Je velmi přesná a stabilní.


### Lotka-Volterra (Model dravec-kořist)
**Princip:**
Soustava dvou nelineárních diferenciálních rovnic prvního řádu.
1. Změna kořisti: Roste rozmnožováním, klesá lovem.
2. Změna dravců: Roste díky jídlu (kořisti), klesá přirozenou smrtí.
Řešení v čase (pomocí Euler nebo RK4) ukazuje fázový posun – nejprve přibude kořist, pak dravci, pak ubyde kořist a nakonec ubydou dravci.