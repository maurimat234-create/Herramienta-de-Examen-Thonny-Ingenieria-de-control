# -*- coding: utf-8 -*-
"""
HERRAMIENTA PARA TECNICOS (version SymPy) - Tarjeta de control de la mesa de corte
Ingenieria de Control - CEUTEC San Pedro Sula - Examen Parcial I

Que hace:
  1. Calcula G_c(s), H_2(s) y H_1 a partir de los componentes de la tarjeta.
  2. Opera funciones de transferencia (serie, paralelo, realimentacion negativa).
  3. Reduce la Figura 1 de forma SIMBOLICA (con Gc, G2, G3, H1, H2) y despues con numeros.
  4. Muestra C(s)/R(s) y Omega(s)/R(s) con coeficientes, polos, ceros y ganancia estacionaria.
  5. Compara contra el resultado hecho a mano (Parte B).
  6. Respuesta al escalon (valor final, sobreimpulso, tiempo de establecimiento).

Requisito: tener sympy instalado.
  En Thonny: menu Herramientas > Administrar paquetes > buscar "sympy" > Instalar.
  (Matplotlib es opcional: solo para la grafica.)

Como usarla: abrir en Thonny y presionar el boton verde (Run / F5).
"""

import cmath
import sympy as sp

s = sp.Symbol("s")

# ======================================================================
# RESULTADO HECHO A MANO (Parte B) -> CAMBIA ESTOS VALORES SI TU d NO ES 8
# ======================================================================
D_MANO = 8
MANO_NUM = 90 * s**2 + 936 * s + 360
MANO_DEN = s**4 + 24 * s**3 + 140 * s**2 + 900 * s + 360


# ======================================================================
# 1. OPERACIONES CON FUNCIONES DE TRANSFERENCIA
#    En sympy una funcion de transferencia es una expresion racional en s.
# ======================================================================
def serie(a, b):
    """SERIE: se multiplican."""
    return a * b


def paralelo(a, b):
    """PARALELO: se suman."""
    return a + b


def realimentar(g, h):
    """LAZO con realimentacion NEGATIVA: G / (1 + G*H)."""
    return g / (1 + g * h)


def como_cociente(expr):
    """Devuelve (numerador, denominador) desarrollados, con denominador monico."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    k = sp.Poly(den, s).LC()                  # coeficiente de la mayor potencia
    return sp.expand(num / k), sp.expand(den / k)


def mostrar(nombre, expr):
    print(f"\n  {nombre} =")
    print(sp.pretty(expr, use_unicode=False))


def mostrar_cociente(nombre, num, den):
    mostrar(nombre, num / den)


def coef(p):
    """Lista de coeficientes de un polinomio en s, mayor potencia primero."""
    return [float(c) for c in sp.Poly(p, s).all_coeffs()]


def fmt_complejo(z):
    z = complex(z)
    if abs(z.imag) < 1e-9:
        return f"{z.real:.6g}"
    signo = "+" if z.imag >= 0 else "-"
    return f"{z.real:.6g} {signo} {abs(z.imag):.6g}j"


# ======================================================================
# 2. ENTRADA DE DATOS CON VALIDACION (numeros exactos con Rational)
# ======================================================================
def pedir(texto, defecto):
    while True:
        r = input(f"  {texto} [{defecto}]: ").strip().replace(",", ".")
        if r == "":
            return sp.Rational(str(defecto))
        try:
            v = sp.Rational(r)
        except (ValueError, TypeError, sp.SympifyError):
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
# 3. RESPUESTA AL ESCALON (por residuos, usando las raices de sympy)
#    c(t) = T(0) + suma de  N(p) / (p * D'(p)) * e^(p t)   (polos simples)
# ======================================================================
def escalon(num, den, polos, K, t_final, pasos=2000):
    dden = sp.diff(den, s)
    N = sp.lambdify(s, num, "math")
    Dp = sp.lambdify(s, dden, "math")
    res = [(complex(p), complex(N(complex(p))) / (complex(p) * complex(Dp(complex(p))))) for p in polos]
    t_list, y_list = [], []
    for k in range(pasos + 1):
        t = t_final * k / pasos
        y = float(K) + sum(r * cmath.exp(p * t) for p, r in res).real
        t_list.append(t)
        y_list.append(y)
    return t_list, y_list


# ======================================================================
# 4. PROGRAMA PRINCIPAL
# ======================================================================
def main():
    titulo("HERRAMIENTA (SymPy) - TARJETA DE CONTROL DE LA MESA DE CORTE")
    print("Presiona Enter para usar el valor que aparece entre [corchetes].")

    d = pedir_digito()
    defectos = {"R1": 10, "R2": 10 + 5 * d, "C": 50, "Ra": 10, "Rf": 10, "Cf": 10,
                "Rin_t": 20, "Rf_t": 10, "Km": 10 + d}

    print("\nControlador PI:")
    R1 = pedir("R1 (kOhm)", defectos["R1"])
    R2 = pedir("R2 (kOhm)", defectos["R2"])
    C = pedir("C  (uF)", defectos["C"])
    print("Filtro del sensor:")
    Ra = pedir("Ra (kOhm)", defectos["Ra"])
    Rf = pedir("Rf (kOhm)", defectos["Rf"])
    Cf = pedir("Cf (uF)", defectos["Cf"])
    print("Tacometro:")
    Rin_t = pedir("R_in (kOhm)", defectos["Rin_t"])
    Rf_t = pedir("R_f  (kOhm)", defectos["Rf_t"])
    print("Motor:")
    Km = pedir("Km", defectos["Km"])

    valores = {"R1": R1, "R2": R2, "C": C, "Ra": Ra, "Rf": Rf, "Cf": Cf,
               "Rin_t": Rin_t, "Rf_t": Rf_t, "Km": Km}
    son_originales = all(valores[k] == sp.Integer(defectos[k]) for k in defectos)

    # ---------------- Parte A: bloques desde los componentes ----------------
    # Unidades: kOhm * uF = milisegundos -> se pasa a segundos con 1e-3
    titulo("1. BLOQUES CALCULADOS DESDE LOS COMPONENTES")
    Kp = R2 / R1                                   # R2/R1 (adimensional)
    Ki = 1 / (R1 * 1000 * C * sp.Rational(1, 10**6))   # 1/(R1*C) con R1 en Ohm y C en F
    Gc_num = Kp + Ki / s
    print(f"  Kp = R2/R1      = {Kp}")
    print(f"  Ki = 1/(R1*C)   = {Ki}  1/s")
    print(f"  Ti = Kp/Ki      = {Kp / Ki}  s")
    mostrar("Gc(s) = Kp + Ki/s", sp.simplify(Gc_num))

    ganancia = Rf / Ra
    tau = Rf * 1000 * Cf * sp.Rational(1, 10**6)   # Rf*Cf en segundos
    H2_num = ganancia / (tau * s + 1)
    print(f"\n  Filtro: ganancia en CC = Rf/Ra = {ganancia} ; constante de tiempo tau = Rf*Cf = {tau} s")
    mostrar("H2(s) = (Rf/Ra)/(Rf*Cf*s + 1)", H2_num)

    h1 = Rf_t / Rin_t
    print(f"\n  H1 = Rf/Rin = {h1}")

    # ---------------- Parte B: reduccion SIMBOLICA ----------------
    titulo("2. REDUCCION DEL DIAGRAMA - FORMA SIMBOLICA")
    Gc, G2, G3, H1, H2 = sp.symbols("Gc G2 G3 H1 H2")
    Gi_sim = realimentar(G2, H1)
    print("  Paso 1 (lazo interno):  Gi = G2 / (1 + G2*H1)")
    mostrar("Gi", Gi_sim)
    G_sim = serie(serie(Gc, Gi_sim), G3)
    print("\n  Paso 2 (serie):  G = Gc * Gi * G3")
    mostrar("G", sp.simplify(G_sim))
    T_sim = sp.simplify(realimentar(G_sim, H2))
    print("\n  Paso 3 (lazo externo):  C/R = G / (1 + G*H2)")
    mostrar("C/R", T_sim)
    W_sim = sp.simplify(realimentar(serie(Gc, Gi_sim), serie(G3, H2)))
    mostrar("Omega/R", W_sim)

    # ---------------- Parte B: con valores ----------------
    titulo("3. REDUCCION CON TUS VALORES")
    G2n = Km / (s + 5)
    G3n = 1 / s
    H1n = h1
    Gi = realimentar(G2n, H1n)
    G = serie(serie(Gc_num, Gi), G3n)
    T = realimentar(G, H2_num)
    W = realimentar(serie(Gc_num, Gi), serie(G3n, H2_num))

    nGi, dGi = como_cociente(Gi)
    mostrar_cociente("Gi(s)", nGi, dGi)
    nT, dT = como_cociente(T)
    mostrar_cociente("C(s)/R(s)", nT, dT)
    nW, dW = como_cociente(W)
    mostrar_cociente("Omega(s)/R(s)", nW, dW)

    # ---------------- Polos, ceros, ganancia ----------------
    titulo("4. POLOS, CEROS Y GANANCIA ESTACIONARIA DE C(s)/R(s)")
    print("  Coeficientes del numerador  :", [round(c, 6) for c in coef(nT)])
    print("  Coeficientes del denominador:", [round(c, 6) for c in coef(dT)])
    polos = sp.Poly(dT, s).nroots()
    ceros = sp.Poly(nT, s).nroots()
    print("  Polos:")
    for z in polos:
        print("     ", fmt_complejo(z))
    print("  Ceros:")
    for z in ceros:
        print("     ", fmt_complejo(z))
    K = sp.simplify(nT.subs(s, 0) / dT.subs(s, 0))
    print(f"  Ganancia en estado estacionario = {K}", "  (OK, es 1)" if K == 1 else "  (revisar)")
    estable = all(complex(z).real < 0 for z in polos)
    print("  Sistema:", "ESTABLE (todos los polos con parte real negativa)" if estable else "INESTABLE")

    # ---------------- Comparacion con el resultado a mano ----------------
    titulo("5. COMPARACION CON EL RESULTADO HECHO A MANO")
    if son_originales and d == D_MANO:
        nM, dM = como_cociente(MANO_NUM / MANO_DEN)
        ok = sp.simplify(nM - nT) == 0 and sp.simplify(dM - dT) == 0
        print("  Numerador")
        print(f"    A mano      : {nM}")
        print(f"    Herramienta : {nT}")
        print("  Denominador")
        print(f"    A mano      : {dM}")
        print(f"    Herramienta : {dT}")
        print("\n  RESULTADO:", "COINCIDE con el calculo a mano" if ok else "NO COINCIDE - revisa tu calculo")
    else:
        print("  (Comparacion omitida: cambiaste algun componente o tu digito no es el de MANO.)")
        print("  Para comparar, usa tus valores originales y edita MANO_NUM / MANO_DEN arriba.")

    # ---------------- Respuesta al escalon ----------------
    titulo("6. RESPUESTA AL ESCALON UNITARIO (posicion del cabezal)")
    if not estable:
        print("  El sistema es inestable; no se simula.")
        return
    lento = min(abs(complex(z).real) for z in polos if abs(complex(z).real) > 1e-9)
    t_final = min(8.0 / lento, 200.0)
    t, y = escalon(nT, dT, polos, K, t_final)
    final = float(K)
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
        print("  (matplotlib no esta instalado: se omite la grafica.)")


if __name__ == "__main__":
    main()
