import numpy as np

# ============================================================
# 1. RELAÇÃO ÁREA-MACH
# ============================================================

def area_mach(M, gamma):
    area_ratio = (1 / M) * (
        (2 / (gamma + 1)) *
        (1 + ((gamma - 1) / 2) * M**2)
    ) ** ((gamma + 1) / (2 * (gamma - 1)))
    return area_ratio


def erro_mach(M, gamma, area_desejada):
    return area_mach(M, gamma) - area_desejada


# ============================================================
# 2. BISSEÇÃO PARA MACH
# ============================================================

def bissecao_mach(gamma, area_desejada, M_min, M_max, tolerancia=1e-6):
    erro_min = erro_mach(M_min, gamma, area_desejada)
    erro_max = erro_mach(M_max, gamma, area_desejada)

    if erro_min * erro_max > 0:
        return None

    while (M_max - M_min) / 2 > tolerancia:
        M_meio = (M_min + M_max) / 2
        erro_meio = erro_mach(M_meio, gamma, area_desejada)

        if erro_meio == 0:
            return M_meio

        if erro_min * erro_meio < 0:
            M_max = M_meio
            erro_max = erro_meio
        else:
            M_min = M_meio
            erro_min = erro_meio

    return (M_min + M_max) / 2


# ============================================================
# 3. RELAÇÕES DE CHOQUE NORMAL
# ============================================================

def mach_pos_choque(M1, gamma):
    M2_squared = (
        (1 + ((gamma - 1) / 2) * M1**2)
        /
        (gamma * M1**2 - ((gamma - 1) / 2))
    )
    return M2_squared**0.5


def pressao_choque(M1, gamma):
    return (1 + (2 * gamma / (gamma + 1)) * (M1**2 - 1))


def densidade_choque(M1, gamma):
    return (((gamma + 1) * M1**2) / ((gamma - 1) * M1**2 + 2))


def temperatura_choque(M1, gamma):
    return pressao_choque(M1, gamma) / densidade_choque(M1, gamma)


def pressao_total_choque(M1, gamma):
    return (
        (((gamma + 1) * M1**2) / ((gamma - 1) * M1**2 + 2)) ** (gamma / (gamma - 1))
        *
        ((gamma + 1) / (2 * gamma * M1**2 - (gamma - 1))) ** (1 / (gamma - 1))
    )


# ============================================================
# 4. GEOMETRIA DO BOCAL
# ============================================================

def area_bocal(x, Ai, At, Ae, Lc, Ld):
    if 0 <= x <= Lc:
        return Ai + (At - Ai) * (x / Lc)
    elif Lc < x <= Lc + Ld:
        return At + (Ae - At) * ((x - Lc) / Ld)
    else:
        return None


# ============================================================
# 5. MACH AO LONGO DO BOCAL
# ============================================================

def mach_bocal(x, Ai, At, Ae, Lc, Ld, gamma):
    A = area_bocal(x, Ai, At, Ae, Lc, Ld)
    if A is None:
        return None

    area_desejada = A / At

    if abs(A - At) < 1e-12:
        return 1.0

    if x < Lc:
        return bissecao_mach(gamma, area_desejada, 0.01, 1.0)
    else:
        return bissecao_mach(gamma, area_desejada, 1.0, 10.0)


# ============================================================
# 6. RELAÇÕES ISENTRÓPICAS
# ============================================================

def pressao_estatica(M, P0, gamma):
    return P0 * (1 + ((gamma - 1) / 2) * M**2) ** (-gamma / (gamma - 1))


def temperatura_estatica(M, T0, gamma):
    return T0 / (1 + ((gamma - 1) / 2) * M**2)


# ============================================================
# 7. ESTADO ATRAVÉS DO CHOQUE
# ============================================================

def estado_choque(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01):
    if x_s is None:
        return None

    As = area_bocal(x_s, Ai, At, Ae, Lc, Ld)
    M1 = mach_bocal(x_s, Ai, At, Ae, Lc, Ld, gamma)

    if As is None or M1 is None or M1 <= 1:
        return None

    P1 = pressao_estatica(M1, P01, gamma)
    T1 = temperatura_estatica(M1, T01, gamma)

    M2 = mach_pos_choque(M1, gamma)
    P2_P1 = pressao_choque(M1, gamma)
    rho2_rho1 = densidade_choque(M1, gamma)
    T2_T1 = temperatura_choque(M1, gamma)

    P2 = P1 * P2_P1
    T2 = T1 * T2_T1

    P02_P01 = pressao_total_choque(M1, gamma)
    P02 = P01 * P02_P01
    T02 = T01

    return {
        "x_s": x_s,
        "A_s": As,
        "M1": M1,
        "P1": P1,
        "T1": T1,
        "M2": M2,
        "P2": P2,
        "T2": T2,
        "P02": P02,
        "T02": T02,
        "P02_P01": P02_P01,
        "P2_P1": P2_P1,
        "rho2_rho1": rho2_rho1,
        "T2_T1": T2_T1
    }


# ============================================================
# 8. ÁREA CRÍTICA PÓS-CHOQUE
# ============================================================

def area_critica_pos_choque(As, M2, gamma):
    area_ratio = area_mach(M2, gamma)
    return As / area_ratio


# ============================================================
# 9. ESTADO NA SAÍDA APÓS O CHOQUE
# ============================================================

def estado_saida_pos_choque(estado, Ae, gamma):
    if estado is None:
        return None

    As = estado["A_s"]
    M2 = estado["M2"]
    P02 = estado["P02"]
    T02 = estado["T02"]

    A2_star = area_critica_pos_choque(As, M2, gamma)
    Ae_A2_star = Ae / A2_star

    Me = bissecao_mach(gamma, Ae_A2_star, 0.01, 1.0)
    if Me is None:
        return None

    Pe = pressao_estatica(Me, P02, gamma)
    Te = temperatura_estatica(Me, T02, gamma)

    return {
        "A2_star": A2_star,
        "Ae_A2_star": Ae_A2_star,
        "Me": Me,
        "Pe": Pe,
        "Te": Te
    }


def pressao_saida_choque(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01):
    estado = estado_choque(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01)
    if estado is None:
        return None
    saida = estado_saida_pos_choque(estado, Ae, gamma)
    if saida is None:
        return None
    return saida["Pe"]


def erro_posicao_choque(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb):
    Pe = pressao_saida_choque(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01)
    if Pe is None:
        return None
    return Pe - Pb


# ============================================================
# 10. LIMITES FÍSICOS DE CONTRAPRESSÃO
# ============================================================

def limites_contrapressao(Ai, At, Ae, Lc, Ld, gamma, P01, T01):
    x_saida = Lc + Ld
    Pb_min = pressao_saida_choque(x_saida, Ai, At, Ae, Lc, Ld, gamma, P01, T01)

    x_garganta = Lc + 1e-6
    Pb_max = pressao_saida_choque(x_garganta, Ai, At, Ae, Lc, Ld, gamma, P01, T01)

    if Pb_min is None or Pb_max is None:
        return None

    return {
        "Pb_min": Pb_min,
        "Pb_max": Pb_max,
        "x_min": x_saida,
        "x_max": x_garganta
    }


def contrapressao_valida(Pb, limites):
    if limites is None:
        return False
    return limites["Pb_min"] <= Pb <= limites["Pb_max"]


# ============================================================
# 11. BISSEÇÃO DA POSIÇÃO DO CHOQUE
# ============================================================

def bissecao_choque(Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb, tolerancia=1e-6):
    limites = limites_contrapressao(Ai, At, Ae, Lc, Ld, gamma, P01, T01)
    if limites is None or not contrapressao_valida(Pb, limites):
        return None

    x_min = Lc + 1e-6
    x_max = Lc + Ld

    erro_min = erro_posicao_choque(x_min, Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb)
    erro_max = erro_posicao_choque(x_max, Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb)

    if erro_min is None or erro_max is None or erro_min * erro_max > 0:
        return None

    while (x_max - x_min) / 2 > tolerancia:
        x_meio = (x_min + x_max) / 2
        erro_meio = erro_posicao_choque(x_meio, Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb)

        if erro_meio is None:
            return None

        if abs(erro_meio) < 1e-12:
            return x_meio

        if erro_min * erro_meio < 0:
            x_max = x_meio
            erro_max = erro_meio
        else:
            x_min = x_meio
            erro_min = erro_meio

    return (x_min + x_max) / 2


# ============================================================
# 12. PERFIL COMPLETO DO BOCAL
# ============================================================

def perfil_bocal(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01, numero_pontos=300):
    if x_s is None:
        return None

    estado = estado_choque(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01)
    if estado is None:
        return None

    A2_star = area_critica_pos_choque(estado["A_s"], estado["M2"], gamma)

    x_antes = np.linspace(0, x_s, numero_pontos // 2, endpoint=False)
    x_depois = np.linspace(x_s, Lc + Ld, numero_pontos // 2)

    x_resultado, A_resultado, M_resultado, P_resultado, T_resultado = [], [], [], [], []

    for x in x_antes:
        A = area_bocal(x, Ai, At, Ae, Lc, Ld)
        M = mach_bocal(x, Ai, At, Ae, Lc, Ld, gamma)
        if M is None: continue
        P = pressao_estatica(M, P01, gamma)
        T = temperatura_estatica(M, T01, gamma)
        x_resultado.append(x); A_resultado.append(A); M_resultado.append(M); P_resultado.append(P); T_resultado.append(T)

    x_resultado.append(x_s); A_resultado.append(estado["A_s"]); M_resultado.append(estado["M1"]); P_resultado.append(estado["P1"]); T_resultado.append(estado["T1"])
    x_resultado.append(x_s); A_resultado.append(estado["A_s"]); M_resultado.append(estado["M2"]); P_resultado.append(estado["P2"]); T_resultado.append(estado["T2"])

    for x in x_depois:
        if abs(x - x_s) < 1e-12: continue
        A = area_bocal(x, Ai, At, Ae, Lc, Ld)
        area_desejada = A / A2_star
        M = bissecao_mach(gamma, area_desejada, 0.01, 1.0)
        if M is None: continue
        P = pressao_estatica(M, estado["P02"], gamma)
        T = temperatura_estatica(M, estado["T02"], gamma)
        x_resultado.append(x); A_resultado.append(A); M_resultado.append(M); P_resultado.append(P); T_resultado.append(T)

    return {
        "x": np.array(x_resultado),
        "A": np.array(A_resultado),
        "M": np.array(M_resultado),
        "P": np.array(P_resultado),
        "T": np.array(T_resultado),
        "A2_star": A2_star,
        "estado_choque": estado
    }


# ============================================================
# 13. FUNÇÃO PRINCIPAL USADA PELA INTERFACE (STREAMLIT)
# ============================================================

def analisar_bocal(Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb):
    limites = limites_contrapressao(Ai, At, Ae, Lc, Ld, gamma, P01, T01)
    if limites is None or not contrapressao_valida(Pb, limites):
        return None

    x_s = bissecao_choque(Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb)
    if x_s is None:
        return None

    perfil = perfil_bocal(x_s, Ai, At, Ae, Lc, Ld, gamma, P01, T01)
    if perfil is None:
        return None

    estado = perfil["estado_choque"]
    saida = estado_saida_pos_choque(estado, Ae, gamma)

    return {
        "xs": x_s,
        "perfil": perfil,
        "estado_choque": estado,
        "saida": saida,
        "limites": limites
    }


# ============================================================
# 14. VALIDAÇÃO AUTOMÁTICA
# ============================================================

def validar_modelo(Pb_min, Pb_max, numero_pontos, Ai, At, Ae, Lc, Ld, gamma, P01, T01):
    resultados = []
    valores_Pb = np.linspace(Pb_min, Pb_max, numero_pontos)

    for Pb in valores_Pb:
        res = analisar_bocal(Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb)
        if res is None:
            resultados.append({
                "Pb": Pb, "xs": None, "As": None, "M1": None, "M2": None,
                "P1": None, "P2": None, "P02": None, "Me": None, "Pe": None,
                "erro_abs": None, "erro_rel": None
            })
            continue

        estado = res["estado_choque"]
        saida = res["saida"]
        Pe = saida["Pe"]
        erro_abs = Pe - Pb
        erro_rel = abs(erro_abs) / Pb * 100

        resultados.append({
            "Pb": Pb,
            "xs": res["xs"],
            "As": estado["A_s"],
            "M1": estado["M1"],
            "M2": estado["M2"],
            "P1": estado["P1"],
            "P2": estado["P2"],
            "P02": estado["P02"],
            "Me": saida["Me"],
            "Pe": Pe,
            "erro_abs": erro_abs,
            "erro_rel": erro_rel
        })

    return resultados
