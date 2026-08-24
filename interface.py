import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from algoritmo import analisar_bocal, validar_modelo, limites_contrapressao

st.set_page_config(page_title="Simulador de Bocal C-D", layout="wide")

st.title("🚀 Simulador de Bocal Convergente-Divergente (C-D)")
st.markdown("Análise de escoamento unidimensional com choque normal interno.")

# ============================================================
# BARRA LATERAL - ENTRADA DE DADOS
# ============================================================
st.sidebar.header("⚙️ Parâmetros de Entrada")

gamma = st.sidebar.number_input("Razão de calores específicos (γ)", value=1.4, step=0.01)
P01 = st.sidebar.number_input("Pressão total na entrada P01 (Pa)", value=101325.0, step=1000.0)
T01 = st.sidebar.number_input("Temperatura total na entrada T01 (K)", value=300.0, step=5.0)

st.sidebar.subheader("Geometria do Bocal")
Ai = st.sidebar.number_input("Área de entrada Ai (m²)", value=2.0, step=0.1)
At = st.sidebar.number_input("Área da garganta At (m²)", value=1.0, step=0.1)
Ae = st.sidebar.number_input("Área de saída Ae (m²)", value=1.6875, step=0.01)

Lc = st.sidebar.number_input("Comprimento convergente Lc (m)", value=1.0, step=0.1)
Ld = st.sidebar.number_input("Comprimento divergente Ld (m)", value=2.0, step=0.1)

# Limites dinâmicos de contrapressão
limites = limites_contrapressao(Ai, At, Ae, Lc, Ld, gamma, P01, T01)

st.sidebar.subheader("Condição de Operação")
if limites:
    pb_min = float(limites["Pb_min"])
    pb_max = float(limites["Pb_max"])
    st.sidebar.info(f"Domínio de choque: {pb_min/1000:.2f} kPa a {pb_max/1000:.2f} kPa")
    Pb = st.sidebar.slider("Contrapressão Pb (Pa)", min_value=pb_min, max_value=pb_max, value=80000.0, step=100.0)
else:
    Pb = st.sidebar.number_input("Contrapressão Pb (Pa)", value=80000.0)

# ============================================================
# PROCESSAMENTO
# ============================================================
resultado = analisar_bocal(Ai, At, Ae, Lc, Ld, gamma, P01, T01, Pb)

if resultado is None:
    st.error("Não foi possível encontrar solução de choque para esses parâmetros ou a contrapressão está fora do domínio físico.")
else:
    perfil = resultado["perfil"]
    estado = resultado["estado_choque"]  # Corrigido de 'estado' para 'estado_choque'
    saida = resultado["saida"]
    xs = resultado["xs"]

    tab1, tab2, tab3 = st.tabs(["📊 Resultados", "📈 Gráficos", "🔍 Validação"])

    # ------------------------------------------------------------
    # TAB 1: RESULTADOS
    # ------------------------------------------------------------
    with tab1:
        st.subheader("Resumo da Posição do Choque")
        col1, col2, col3 = st.columns(3)
        col1.metric("Posição do Choque (xs)", f"{xs:.4f} m")
        col2.metric("Área no Choque (As)", f"{estado['A_s']:.4f} m²")
        col3.metric("Mach Antes do Choque (M1)", f"{estado['M1']:.3f}")

        st.markdown("---")
        c1, c2 = st.columns(2)

        with c1:
            st.subheader("Propriedades Pré e Pós Choque")
            df_choque = pd.DataFrame({
                "Propriedade": ["Mach (M)", "Pressão Estática (Pa)", "Temperatura (K)"],
                "Antes (1)": [f"{estado['M1']:.4f}", f"{estado['P1']:.2f}", f"{estado['T1']:.2f}"],
                "Depois (2)": [f"{estado['M2']:.4f}", f"{estado['P2']:.2f}", f"{estado['T2']:.2f}"]
            })
            st.table(df_choque)

        with c2:
            st.subheader("Condições na Saída do Bocal")
            df_saida = pd.DataFrame({
                "Propriedade": ["Mach na Saída (Me)", "Pressão de Saída (Pe)", "Temperatura de Saída (Te)"],
                "Valor": [f"{saida['Me']:.4f}", f"{saida['Pe']:.2f} Pa", f"{saida['Te']:.2f} K"]
            })
            st.table(df_saida)

    # ------------------------------------------------------------
    # TAB 2: GRÁFICOS
    # ------------------------------------------------------------
    with tab2:
        st.subheader("Distribuição ao longo do Bocal")

        fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

        # Gráfico de Mach
        ax[0].plot(perfil["x"], perfil["M"], label="Mach", color="blue")
        ax[0].axvline(xs, color="red", linestyle="--", label="Onda de Choque")
        ax[0].set_ylabel("Número de Mach")
        ax[0].grid(True)
        ax[0].legend()

        # Gráfico de Pressão
        ax[1].plot(perfil["x"], perfil["P"] / 1000, label="Pressão Estática", color="green")
        ax[1].axhline(Pb / 1000, color="orange", linestyle=":", label="Contrapressão Pb")
        ax[1].axvline(xs, color="red", linestyle="--")
        ax[1].set_ylabel("Pressão (kPa)")
        ax[1].grid(True)
        ax[1].legend()

        # Gráfico de Temperatura
        ax[2].plot(perfil["x"], perfil["T"], label="Temperatura Estática", color="red")
        ax[2].axvline(xs, color="red", linestyle="--")
        ax[2].set_xlabel("Posição x (m)")
        ax[2].set_ylabel("Temperatura (K)")
        ax[2].grid(True)
        ax[2].legend()

        st.pyplot(fig)

    # ------------------------------------------------------------
    # TAB 3: VALIDAÇÃO
    # ------------------------------------------------------------
    with tab3:
        st.subheader("Tabela de Validação de Contrapressões")
        res_validacao = validar_modelo(pb_min, pb_max, 10, Ai, At, Ae, Lc, Ld, gamma, P01, T01)
        df_val = pd.DataFrame(res_validacao)
        st.dataframe(df_val)
