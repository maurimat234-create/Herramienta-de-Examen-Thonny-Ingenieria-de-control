# -*- coding: utf-8 -*-
"""
HERRAMIENTA PARA TECNICOS - Tarjeta de control de la mesa de corte
Ingenieria de Control - CEUTEC San Pedro Sula - Examen Parcial I

Que hace:
  1. Calcula G_c(s), H_2(s) y H_1 a partir de los componentes de la tarjeta.
  2. Opera funciones de transferencia (serie, paralelo, realimentacion negativa).
  3. Reduce el diagrama de la Figura 1 y muestra C(s)/R(s) y Omega(s)/R(s)
     con coeficientes, polos, ceros y ganancia en estado estacionario.
  4. Compara contra el resultado hecho a mano (Parte B).
  5. Calcula la respuesta al escalon (valor final, sobreimpulso, tiempo de establecimiento).

Como usarla: abrir en Thonny y presionar el boton verde (Run / F5).
No necesita instalar nada. Si tienes matplotlib, ademas dibuja la grafica.
"""

# ======================================================================
# RESULTADO HECHO A MANO (Parte B) -> CAMBIA ESTOS VALORES SI TU d NO ES 8
# Coeficientes de mayor a menor potencia de s.
# ======================================================================
D_MANO = 8
MANO_NUM = [90, 936, 360]            # numerador de C(s)/R(s) a mano
MANO_DEN = [1, 24, 140, 900, 360]    # denominador de C(s)/R(s) a mano


# ======================================================================
# 1. POLINOMIOS: lista de coeficientes, potencia mayor primero
#    Ejemplo: 3s + 2  ->  [3, 2]
# ======================================================================
def limpiar(p, tol=1e-12):
    """Quita ceros sobrantes a la izquierda."""
    p = list(p)
    while len(p) > 1 and abs(p[0]) < tol:
        p.pop(0)
    return p


def p_suma(a, b):
    """Suma dos polinomios (alinea por la potencia mas baja)."""
    n = max(len(a), len(b))
    a = [0.0] * (n - len(a)) + list(a)
    b = [0.0] * (n - len(b)) + list(b)
    return limpiar([x + y for x, y in zip(a, b)])


def p_mul(a, b):
    """Multiplica dos polinomios."""
    r = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return limpiar(r)


def p_eval(p, x):
    """Evalua el polinomio en x (metodo de Horner)."""
    r = 0
    for c in p:
        r = r * x + c
    return r


def poli_str(p, var="s"):
    """Convierte [90, 936, 360] en '90 s^2 + 936 s + 360'."""
    n = len(p) - 1
    texto = ""
    for i, c in enumerate(p):
        e = n - i
        if abs(c) < 1e-12:
            continue
        mag = f"{abs(c):.6g}"
        if e == 0:
            term = mag
        else:
            coef = "" if abs(abs(c) - 1) < 1e-12 else mag + " "
            term = f"{coef}{var}" + (f"^{e}" if e > 1 else "")
        if texto == "":
            texto = ("-" if c < 0 else "") + term
        else:
            texto += (" - " if c < 0 else " + ") + term
    return texto if texto else "0"


def raices(p):
    """Raices de un polinomio. Usa numpy si existe; si no, Durand-Kerner."""
    p = limpiar(p)
    n = len(p) - 1
    if n < 1:
        return []
    try:
        import numpy as np
        return [complex(z) for z in np.roots(p)]
    except ImportError:
        a = [c / p[0] for c in p]
        r = [(0.4 + 0.9j) ** k for k in range(n)]      # estimados iniciales
        for _ in range(1000):
            for i in range(n):
                val = sum(a[k] * r[i] ** (n - k) for k in range(n + 1))
                den = 1
                for j in range(n):
                    if j != i:
                        den *= (r[i] - r[j])
                r[i] -= val / den
        return r


def fmt_complejo(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:.6g}"
    signo = "+" if z.imag >= 0 else "-"
    return f"{z.real:.6g} {signo} {abs(z.imag):.6g}j"


# ======================================================================
# 2. FUNCION DE TRANSFERENCIA = numerador / denominador
# ======================================================================
class TF:
    def __init__(self, num, den):
        self.num = limpiar(num)
        self.den = limpiar(den)

    def __mul__(self, otro):
        """SERIE: se multiplican numeradores y denominadores."""
        return TF(p_mul(self.num, otro.num), p_mul(self.den, otro.den))

    def __add__(self, otro):
        """PARALELO: a/b + c/d = (a*d + c*b) / (b*d)."""
        num = p_suma(p_mul(self.num, otro.den), p_mul(otro.num, self.den))
        return TF(num, p_mul(self.den, otro.den))

    def realimentar(self, H):
        """LAZO con realimentacion NEGATIVA: G / (1 + G*H)."""
        num = p_mul(self.num, H.den)
        den = p_suma(p_mul(self.den, H.den), p_mul(self.num, H.num))
        return TF(num, den)

    def normalizar(self):
        """Divide todo entre el coeficiente mayor del denominador (queda monico)."""
        k = self.den[0]
        return TF([c / k for c in self.num], [c / k for c in self.den])

    def ganancia_cc(self):
        """Ganancia en estado estacionario: valor de G(s) cuando s = 0."""
        d0 = p_eval(self.den, 0)
        if abs(d0) < 1e-12:
            return float("inf")
        return p_eval(self.num, 0) / d0

    def mostrar(self, nombre):
        n, d = poli_str(self.num), poli_str(self.den)
        ancho = max(len(n), len(d)) + 2
        print(f"  {nombre} =")
        print("      " + n.center(ancho))
        print("      " + "-" * ancho)
        print("      " + d.center(ancho))


# ======================================================================
# 3. BLOQUES DESDE LOS COMPONENTES (Parte A)
# ======================================================================
def bloque_PI(R1, R2, C):
    """Etapa PI + inversor:  G_c = R2/R1 + 1/(R1*C) * 1/s"""
    Kp = R2 / R1
    Ki = 1.0 / (R1 * C)
    return Kp, Ki, TF([Kp, Ki], [1, 0])        # (Kp*s + Ki) / s


def bloque_filtro(Ra, Rf, Cf):
    """Filtro + inversor:  H_2 = (Rf/Ra) / (Rf*Cf*s + 1)"""
    ganancia = Rf / Ra
    tau = Rf * Cf
    return ganancia, tau, TF([ganancia], [tau, 1])


def bloque_tacometro(Rin, Rf):
    """Inversor + inversor:  H_1 = Rf / Rin"""
    return Rf / Rin


# ======================================================================
# 4. RESPUESTA AL ESCALON (simulacion numerica, Runge-Kutta 4)
# ======================================================================
def escalon(T, t_final, pasos=4000):
    """Devuelve (tiempos, salida) para entrada escalon unitario."""
    T = T.normalizar()
    n = len(T.den) - 1
    a = T.den[1:]                                   # a1..an
    b = [0.0] * (n + 1 - len(T.num)) + T.num        # b0..bn
    d0 = b[0]
    c = [b[k] - a[k - 1] * d0 for k in range(1, n + 1)]
    C = c[::-1]                                     # forma canonica controlable

    def deriv(x):
        # x' = A x + B u, con u = 1 (escalon unitario)
        ultimo = -sum(a[n - 1 - i] * x[i] for i in range(n)) + 1.0
        return x[1:] + [ultimo]

    x = [0.0] * n
    dt = t_final / pasos
    t_list, y_list = [], []
    for k in range(pasos + 1):
        y_list.append(sum(C[i] * x[i] for i in range(n)) + d0)
        t_list.append(k * dt)
        k1 = deriv(x)
        k2 = deriv([x[i] + dt / 2 * k1[i] for i in range(n)])
        k3 = deriv([x[i] + dt / 2 * k2[i] for i in range(n)])
        k4 = deriv([x[i] + dt * k3[i] for i in range(n)])
        x = [x[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(n)]
    return t_list, y_list


# ======================================================================
# 5. ENTRADA DE DATOS CON VALIDACION
# ======================================================================
def pedir(texto, defecto):
    while True:
        r = input(f"  {texto} [{defecto:g}]: ").strip().replace(",", ".")
        if r == "":
            return float(defecto)
        try:
            v = float(r)
        except ValueError:
            print("    -> Escribe un numero (ejemplo: 50 o 2.5).")
            continue
        if v <= 0:
            print("    -> Debe ser mayor que 0.")
            continue
        return v


def pedir_digito():
    while True:
        r = input("  Ultimo digito de tu numero de cuenta (0-9): ").strip()
        if r.isdigit() and len(r) == 1:
            return int(r)
        print("    -> Escribe un solo digito entre 0 y 9.")


def titulo(texto):
    print("\n" + "=" * 62)
    print(texto)
    print("=" * 62)


# ======================================================================
# 6. PROGRAMA PRINCIPAL
# ======================================================================
def main():
    titulo("HERRAMIENTA - TARJETA DE CONTROL DE LA MESA DE CORTE")
    print("Presiona Enter para usar el valor que aparece entre [corchetes].")

    d = pedir_digito()
    defectos = {
        "R1": 10.0, "R2": 10.0 + 5 * d, "C": 50.0,
        "Ra": 10.0, "Rf": 10.0, "Cf": 10.0,
        "Rin_t": 20.0, "Rf_t": 10.0, "Km": 10.0 + d,
    }

    print("\nControlador PI:")
    R1 = pedir("R1 (kOhm)", defectos["R1"]) * 1e3
    R2 = pedir("R2 (kOhm)", defectos["R2"]) * 1e3
    C = pedir("C  (uF)", defectos["C"]) * 1e-6
    print("Filtro del sensor:")
    Ra = pedir("Ra (kOhm)", defectos["Ra"]) * 1e3
    Rf = pedir("Rf (kOhm)", defectos["Rf"]) * 1e3
    Cf = pedir("Cf (uF)", defectos["Cf"]) * 1e-6
    print("Tacometro:")
    Rin_t = pedir("R_in (kOhm)", defectos["Rin_t"]) * 1e3
    Rf_t = pedir("R_f  (kOhm)", defectos["Rf_t"]) * 1e3
    print("Motor:")
    Km = pedir("Km", defectos["Km"])

    valores = {"R1": R1 / 1e3, "R2": R2 / 1e3, "C": C / 1e-6, "Ra": Ra / 1e3,
               "Rf": Rf / 1e3, "Cf": Cf / 1e-6, "Rin_t": Rin_t / 1e3,
               "Rf_t": Rf_t / 1e3, "Km": Km}
    son_originales = all(abs(valores[k] - defectos[k]) < 1e-9 for k in defectos)

    # ---------------- Parte A: bloques desde los componentes ----------------
    titulo("1. BLOQUES CALCULADOS DESDE LOS COMPONENTES")
    Kp, Ki, Gc = bloque_PI(R1, R2, C)
    print(f"  Kp = R2/R1      = {Kp:.6g}")
    print(f"  Ki = 1/(R1*C)   = {Ki:.6g}  1/s")
    print(f"  Ti = Kp/Ki      = {Kp / Ki:.6g}  s")
    Gc.mostrar("Gc(s)")
    gan, tau, H2 = bloque_filtro(Ra, Rf, Cf)
    print(f"\n  Filtro: ganancia en CC = Rf/Ra = {gan:.6g} ; constante de tiempo tau = Rf*Cf = {tau:.6g} s")
    H2.mostrar("H2(s)")
    h1 = bloque_tacometro(Rin_t, Rf_t)
    print(f"\n  H1 = Rf/Rin = {h1:.6g}")

    # ---------------- Parte B: reduccion del diagrama ----------------
    titulo("2. REDUCCION DEL DIAGRAMA DE LA FIGURA 1")
    G2 = TF([Km], [1, 5])           # motor
    G3 = TF([1], [1, 0])            # integrador (velocidad -> posicion)
    H1 = TF([h1], [1])              # tacometro como bloque constante

    Gi = G2.realimentar(H1)         # paso 1: lazo interno
    print("  Paso 1 (lazo interno): Gi = G2 / (1 + G2*H1)")
    Gi.normalizar().mostrar("Gi(s)")

    G = Gc * Gi * G3                # paso 2: serie
    print("\n  Paso 2 (serie): G = Gc * Gi * G3")
    G.normalizar().mostrar("G(s)")

    T = G.realimentar(H2).normalizar()          # paso 3: lazo externo
    print("\n  Paso 3 (lazo externo): C/R = G / (1 + G*H2)")
    T.mostrar("C(s)/R(s)")

    # Omega/R: camino directo Gc*Gi ; realimentacion G3*H2
    W = (Gc * Gi).realimentar(G3 * H2).normalizar()
    print("\n  Velocidad del motor: Omega/R = Gc*Gi / (1 + Gc*Gi*G3*H2)")
    W.mostrar("Omega(s)/R(s)")

    # ---------------- Polos, ceros, ganancia ----------------
    titulo("3. POLOS, CEROS Y GANANCIA ESTACIONARIA DE C(s)/R(s)")
    print("  Coeficientes del numerador  :", [round(c, 6) for c in T.num])
    print("  Coeficientes del denominador:", [round(c, 6) for c in T.den])
    polos = raices(T.den)
    ceros = raices(T.num)
    print("  Polos:")
    for z in polos:
        print("     ", fmt_complejo(z))
    print("  Ceros:")
    for z in ceros:
        print("     ", fmt_complejo(z))
    K = T.ganancia_cc()
    print(f"  Ganancia en estado estacionario = {K:.6g}", "  (OK, es 1)" if abs(K - 1) < 1e-9 else "  (revisar)")
    estable = all(z.real < 0 for z in polos)
    print("  Sistema:", "ESTABLE (todos los polos con parte real negativa)" if estable else "INESTABLE")

    # ---------------- Comparacion con el resultado a mano ----------------
    titulo("4. COMPARACION CON EL RESULTADO HECHO A MANO")
    if son_originales and d == D_MANO:
        # Se normaliza el resultado a mano igual que el de la herramienta
        mano = TF(MANO_NUM, MANO_DEN).normalizar()
        ok_n = len(mano.num) == len(T.num) and all(abs(x - y) < 1e-6 for x, y in zip(mano.num, T.num))
        ok_d = len(mano.den) == len(T.den) and all(abs(x - y) < 1e-6 for x, y in zip(mano.den, T.den))
        print(f"  {'':<4}{'A mano':>20}{'Herramienta':>16}")
        print("  Numerador:")
        for a_, b_ in zip(mano.num, T.num):
            print(f"    {a_:>24.6f}{b_:>16.6f}")
        print("  Denominador:")
        for a_, b_ in zip(mano.den, T.den):
            print(f"    {a_:>24.6f}{b_:>16.6f}")
        print("\n  RESULTADO:", "COINCIDE con el calculo a mano" if (ok_n and ok_d) else "NO COINCIDE - revisa tu calculo")
    else:
        print("  (Comparacion omitida: cambiaste algun componente o tu digito no es el de MANO.)")
        print("  Para comparar, usa tus valores originales y edita MANO_NUM / MANO_DEN arriba.")

    # ---------------- Respuesta al escalon ----------------
    titulo("5. RESPUESTA AL ESCALON UNITARIO (posicion del cabezal)")
    if not estable:
        print("  El sistema es inestable; no se simula.")
        return
    lento = min(abs(z.real) for z in polos if abs(z.real) > 1e-9)
    t_final = min(8.0 / lento, 200.0)
    t, y = escalon(T, t_final)
    final = K
    sobre = max(0.0, (max(y) - final) / abs(final) * 100)
    ts = 0.0
    for ti, yi in zip(t, y):
        if abs(yi - final) > 0.02 * abs(final):
            ts = ti
    print(f"  Valor final             : {y[-1]:.4f}   (esperado {final:.4f})")
    print(f"  Sobreimpulso maximo     : {sobre:.2f} %")
    print(f"  Tiempo de establecimiento (2 %): {ts:.3f} s")
    print("  Muestras (t, c(t)):")
    pasos = len(t) - 1
    for k_ in range(0, pasos + 1, pasos // 10):
        print(f"     t = {t[k_]:8.3f} s   c = {y[k_]:.4f}")
    try:
        import matplotlib.pyplot as plt
        plt.plot(t, y, label="c(t)")
        plt.axhline(final, linestyle="--", color="gray", label="valor final")
        plt.title("Respuesta al escalon de C(s)/R(s)")
        plt.xlabel("t (s)")
        plt.ylabel("posicion")
        plt.grid(True)
        plt.legend()
        plt.show()
    except ImportError:
        print("  (matplotlib no esta instalado: se omite la grafica. En Thonny: Herramientas > Administrar paquetes.)")


if __name__ == "__main__":
    main()
    