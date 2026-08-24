import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from algoritmo import analisar_bocal, validar_modelo


# ============================================================
# CONFIGURAÇÃO DA PÁGINA
# ============================================================

st.set_page_config(
    page_title="Analisador de Bocais C-D",
    page_icon="🚀",
    layout="wide"
)


# ============================================================
# CSS
# ============================================================

st.markdown(
    """
    <style>

    .main {
        padding-top: 1rem;
    }

    .titulo-principal {
        font-size: 2.4rem;
        font-weight: 700;
        margin-bottom: 0rem;
    }

    .subtitulo {
        font-size: 1.05rem;
        color: #777;
        margin-bottom: 1.5rem;
    }

    .card {
        padding: 1rem;
        border-radius: 10px;
        border: 1px solid #ddd;
        background-color: #fafafa;
        margin-bottom: 1rem;
    }

    .resultado-destaque {
        font-size: 1.8rem;
        font-weight: 700;
    }

    </style>
    """,
    unsafe_allow_html=True
)


# ============================================================
# CABEÇALHO
# ============================================================

st.markdown(
    '<div class="titulo-principal">'
    'Analisador de Bocais Convergente-Divergentes'
    '</div>',
    unsafe_allow_html=True
)

st.markdown(
    '<div class="subtitulo">'
    'Modelo analítico unidimensional para escoamento compressível '
    'com onda de choque normal'
    '</div>',
    unsafe_allow_html=True
)


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("Parâmetros do modelo")


# ------------------------------------------------------------
# Escoamento
# ------------------------------------------------------------

st.sidebar.subheader("Condições de estagnação")

gamma = st.sidebar.number_input(
    "γ — Razão de calores específicos",
    min_value=1.01,
    max_value=2.00,
    value=1.40,
    step=0.01,
    format="%.2f"
)

P01 = st.sidebar.number_input(
    "P₀₁ — Pressão total (Pa)",
    min_value=1.0,
    value=101325.0,
    step=1000.0,
    format="%.1f"
)

T01 = st.sidebar.number_input(
    "T₀₁ — Temperatura total (K)",
    min_value=1.0,
    value=300.0,
    step=10.0,
    format="%.1f"
)


# ------------------------------------------------------------
# Geometria
# ------------------------------------------------------------

st.sidebar.subheader("Geometria do bocal")

Ai = st.sidebar.number_input(
    "Aᵢ — Área de entrada (m²)",
    min_value=0.001,
    value=2.000,
    step=0.1,
    format="%.4f"
)

At = st.sidebar.number_input(
    "Aₜ — Área da garganta (m²)",
    min_value=0.001,
    value=1.000,
    step=0.1,
    format="%.4f"
)

Ae = st.sidebar.number_input(
    "Aₑ — Área de saída (m²)",
    min_value=0.001,
    value=1.6875,
    step=0.1,
    format="%.4f"
)

Lc = st.sidebar.number_input(
    "Comprimento convergente (m)",
    min_value=0.001,
    value=1.000,
    step=0.1,
    format="%.4f"
)

Ld = st.sidebar.number_input(
    "Comprimento divergente (m)",
    min_value=0.001,
    value=2.000,
    step=0.1,
    format="%.4f"
)


# ------------------------------------------------------------
# Contrapressão
# ------------------------------------------------------------

st.sidebar.subheader("Condição de saída")

Pb = st.sidebar.number_input(
    "Pᵦ — Contrapressão (Pa)",
    min_value=1.0,
    value=75000.0,
    step=1000.0,
    format="%.1f"
)


# ============================================================
# BOTÃO DE CÁLCULO
# ============================================================

calcular = st.sidebar.button(
    "CALCULAR CASO",
    use_container_width=True
)


# ============================================================
# VALIDAÇÕES BÁSICAS
# ============================================================

if Ai <= At:

    st.sidebar.warning(
        "A área de entrada deve ser maior que a garganta."
    )

if Ae <= At:

    st.sidebar.warning(
        "A área de saída deve ser maior que a garganta."
    )


# ============================================================
# ESTADO DA APLICAÇÃO
# ============================================================

if "resultado" not in st.session_state:

    st.session_state.resultado = None


# ============================================================
# CÁLCULO
# ============================================================

if calcular:

    if Ai <= At:

        st.error(
            "A área de entrada deve ser maior que a área da garganta."
        )

        st.stop()

    if Ae <= At:

        st.error(
            "A área de saída deve ser maior que a área da garganta."
        )

        st.stop()

    resultado = analisar_bocal(
        gamma=gamma,
        P01=P01,
        T01=T01,
        Ai=Ai,
        At=At,
        Ae=Ae,
        Lc=Lc,
        Ld=Ld,
        Pb=Pb
    )

    st.session_state.resultado = resultado


# ============================================================
# RECUPERA RESULTADO
# ============================================================

resultado = st.session_state.resultado


# ============================================================
# SE NÃO HOUVER RESULTADO
# ============================================================

if resultado is None:

    st.info(
        "Configure os parâmetros do bocal no painel lateral "
        "e clique em **CALCULAR CASO**."
    )

    st.markdown("### Modelo")

    col1, col2, col3 = st.columns(3)

    with col1:

        st.metric(
            "Modelo",
            "Quase-1D"
        )

    with col2:

        st.metric(
            "Choque",
            "Normal"
        )

    with col3:

        st.metric(
            "Hipóteses",
            "Adiabático / invíscido"
        )

    st.stop()


# ============================================================
# EXTRAI RESULTADOS
# ============================================================

estado = resultado["estado"]

saida = resultado["saida"]

perfil = resultado["perfil"]

x_s = resultado["x_s"]

erro_abs = resultado["erro_abs"]

erro_rel = resultado["erro_rel"]


# ============================================================
# INDICADORES PRINCIPAIS
# ============================================================

st.success("Cálculo concluído com sucesso.")


st.subheader("Resultado principal")


col1, col2, col3, col4 = st.columns(4)


with col1:

    st.metric(
        "Posição do choque",
        f"{x_s:.6f} m"
    )


with col2:

    st.metric(
        "Área no choque",
        f"{estado['A_s']:.6f} m²"
    )


with col3:

    st.metric(
        "Mach antes do choque",
        f"{estado['M1']:.6f}"
    )


with col4:

    st.metric(
        "Mach depois do choque",
        f"{estado['M2']:.6f}"
    )


# ============================================================
# VERIFICAÇÃO Pe = Pb
# ============================================================

st.subheader("Verificação da condição de saída")


col1, col2, col3 = st.columns(3)


with col1:

    st.metric(
        "Contrapressão Pb",
        f"{Pb / 1000:.6f} kPa"
    )


with col2:

    st.metric(
        "Pressão calculada Pe",
        f"{saida['Pe'] / 1000:.6f} kPa"
    )


with col3:

    st.metric(
        "Erro Pe − Pb",
        f"{erro_abs:.6f} Pa"
    )


st.caption(
    f"Erro relativo: {erro_rel:.8f} %"
)


# ============================================================
# ABAS
# ============================================================

aba_resultados, aba_graficos, aba_validacao = st.tabs(
    [
        "Resultados",
        "Gráficos",
        "Validação"
    ]
)


# ============================================================
# ABA RESULTADOS
# ============================================================

with aba_resultados:

    st.subheader("Estado antes do choque")

    col1, col2, col3, col4 = st.columns(4)

    with col1:

        st.metric(
            "M₁",
            f"{estado['M1']:.6f}"
        )

    with col2:

        st.metric(
            "P₁",
            f"{estado['P1'] / 1000:.6f} kPa"
        )

    with col3:

        st.metric(
            "T₁",
            f"{estado['T1']:.6f} K"
        )

    with col4:

        st.metric(
            "Aₛ",
            f"{estado['A_s']:.6f} m²"
        )


    st.subheader("Estado depois do choque")

    col1, col2, col3 = st.columns(3)

    with col1:

        st.metric(
            "M₂",
            f"{estado['M2']:.6f}"
        )

    with col2:

        st.metric(
            "P₂",
            f"{estado['P2'] / 1000:.6f} kPa"
        )

    with col3:

        st.metric(
            "T₂",
            f"{estado['T2']:.6f} K"
        )


    st.subheader("Condições totais")

    col1, col2, col3, col4 = st.columns(4)

    with col1:

        st.metric(
            "P₀₁",
            f"{P01 / 1000:.3f} kPa"
        )

    with col2:

        st.metric(
            "P₀₂",
            f"{estado['P02'] / 1000:.3f} kPa"
        )

    with col3:

        st.metric(
            "P₀₂ / P₀₁",
            f"{estado['P02_P01']:.6f}"
        )

    with col4:

        st.metric(
            "T₀₂",
            f"{estado['T02']:.3f} K"
        )


    st.subheader("Estado na saída")

    col1, col2, col3, col4 = st.columns(4)

    with col1:

        st.metric(
            "A₂*",
            f"{saida['A2_star']:.6f} m²"
        )

    with col2:

        st.metric(
            "Aₑ / A₂*",
            f"{saida['Ae_A2_star']:.6f}"
        )

    with col3:

        st.metric(
            "Mₑ",
            f"{saida['Me']:.6f}"
        )

    with col4:

        st.metric(
            "Tₑ",
            f"{saida['Te']:.6f} K"
        )


    # --------------------------------------------------------
    # Tabela resumida
    # --------------------------------------------------------

    st.subheader("Resumo")

    dados = {

        "Grandeza": [
            "Contrapressão",
            "Posição do choque",
            "Área do choque",
            "Mach antes do choque",
            "Pressão antes do choque",
            "Temperatura antes do choque",
            "Mach depois do choque",
            "Pressão depois do choque",
            "Temperatura depois do choque",
            "Pressão total após choque",
            "Mach na saída",
            "Pressão na saída",
            "Temperatura na saída",
            "Erro absoluto",
            "Erro relativo"
        ],

        "Valor": [

            f"{Pb:.6f} Pa",

            f"{x_s:.6f} m",

            f"{estado['A_s']:.6f} m²",

            f"{estado['M1']:.6f}",

            f"{estado['P1']:.6f} Pa",

            f"{estado['T1']:.6f} K",

            f"{estado['M2']:.6f}",

            f"{estado['P2']:.6f} Pa",

            f"{estado['T2']:.6f} K",

            f"{estado['P02']:.6f} Pa",

            f"{saida['Me']:.6f}",

            f"{saida['Pe']:.6f} Pa",

            f"{saida['Te']:.6f} K",

            f"{erro_abs:.6f} Pa",

            f"{erro_rel:.8f} %"
        ]
    }

    tabela_resumo = pd.DataFrame(dados)

    st.dataframe(
        tabela_resumo,
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# ABA GRÁFICOS
# ============================================================

with aba_graficos:

    st.subheader("Distribuição ao longo do bocal")


    # --------------------------------------------------------
    # Gráfico de Mach
    # --------------------------------------------------------

    fig_mach, ax_mach = plt.subplots(
        figsize=(10, 5)
    )

    ax_mach.plot(
        perfil["x"],
        perfil["M"],
        linewidth=2
    )

    ax_mach.axvline(
        x_s,
        linestyle="--",
        linewidth=1.5,
        label="Onda de choque"
    )

    ax_mach.set_xlabel(
        "Posição x (m)"
    )

    ax_mach.set_ylabel(
        "Número de Mach"
    )

    ax_mach.set_title(
        "Distribuição do número de Mach"
    )

    ax_mach.grid(
        alpha=0.3
    )

    ax_mach.legend()

    st.pyplot(
        fig_mach,
        use_container_width=True
    )

    plt.close(fig_mach)


    # --------------------------------------------------------
    # Gráfico de pressão
    # --------------------------------------------------------

    fig_pressao, ax_pressao = plt.subplots(
        figsize=(10, 5)
    )

    ax_pressao.plot(
        perfil["x"],
        perfil["P"] / 1000,
        linewidth=2,
        label="Pressão estática"
    )

    ax_pressao.axvline(
        x_s,
        linestyle="--",
        linewidth=1.5,
        label="Onda de choque"
    )

    ax_pressao.axhline(
        Pb / 1000,
        linestyle=":",
        linewidth=1.5,
        label="Contrapressão"
    )

    ax_pressao.set_xlabel(
        "Posição x (m)"
    )

    ax_pressao.set_ylabel(
        "Pressão (kPa)"
    )

    ax_pressao.set_title(
        "Distribuição da pressão estática"
    )

    ax_pressao.grid(
        alpha=0.3
    )

    ax_pressao.legend()

    st.pyplot(
        fig_pressao,
        use_container_width=True
    )

    plt.close(fig_pressao)


    # --------------------------------------------------------
    # Gráfico de temperatura
    # --------------------------------------------------------

    fig_temp, ax_temp = plt.subplots(
        figsize=(10, 5)
    )

    ax_temp.plot(
        perfil["x"],
        perfil["T"],
        linewidth=2
    )

    ax_temp.axvline(
        x_s,
        linestyle="--",
        linewidth=1.5,
        label="Onda de choque"
    )

    ax_temp.set_xlabel(
        "Posição x (m)"
    )

    ax_temp.set_ylabel(
        "Temperatura (K)"
    )

    ax_temp.set_title(
        "Distribuição da temperatura estática"
    )

    ax_temp.grid(
        alpha=0.3
    )

    ax_temp.legend()

    st.pyplot(
        fig_temp,
        use_container_width=True
    )

    plt.close(fig_temp)


    # --------------------------------------------------------
    # Geometria do bocal
    # --------------------------------------------------------

    st.subheader("Geometria do bocal")

    x_geom = np.linspace(
        0,
        Lc + Ld,
        300
    )

    A_geom = np.array(
        [
            (
                Ai + (At - Ai) * (x / Lc)
                if x <= Lc
                else
                At + (Ae - At) * ((x - Lc) / Ld)
            )
            for x in x_geom
        ]
    )

    fig_geom, ax_geom = plt.subplots(
        figsize=(10, 4)
    )

    ax_geom.plot(
        x_geom,
        A_geom,
        linewidth=2
    )

    ax_geom.axvline(
        x_s,
        linestyle="--",
        linewidth=1.5,
        label="Onda de choque"
    )

    ax_geom.set_xlabel(
        "Posição x (m)"
    )

    ax_geom.set_ylabel(
        "Área (m²)"
    )

    ax_geom.set_title(
        "Perfil de área do bocal"
    )

    ax_geom.grid(
        alpha=0.3
    )

    ax_geom.legend()

    st.pyplot(
        fig_geom,
        use_container_width=True
    )

    plt.close(fig_geom)


# ============================================================
# ABA VALIDAÇÃO
# ============================================================

with aba_validacao:

    st.subheader(
        "Varredura da contrapressão"
    )

    st.write(
        "Calcule automaticamente a posição da onda de choque "
        "para uma faixa de contrapressões."
    )


    col1, col2, col3 = st.columns(3)


    with col1:

        Pb_min = st.number_input(
            "Pb mínimo (Pa)",
            min_value=1.0,
            value=60000.0,
            step=1000.0
        )


    with col2:

        Pb_max = st.number_input(
            "Pb máximo (Pa)",
            min_value=1.0,
            value=90000.0,
            step=1000.0
        )


    with col3:

        numero_pontos = st.number_input(
            "Número de pontos",
            min_value=2,
            max_value=1000,
            value=13,
            step=1
        )


    executar_validacao = st.button(
        "EXECUTAR VARREDURA",
        use_container_width=True
    )


    if executar_validacao:

        if Pb_min >= Pb_max:

            st.error(
                "Pb mínimo deve ser menor que Pb máximo."
            )

        else:

            with st.spinner(
                "Executando a varredura..."
            ):

                resultados = validar_modelo(
                    Pb_min=Pb_min,
                    Pb_max=Pb_max,
                    numero_pontos=int(numero_pontos),
                    Ai=Ai,
                    At=At,
                    Ae=Ae,
                    Lc=Lc,
                    Ld=Ld,
                    gamma=gamma,
                    P01=P01,
                    T01=T01
                )


            # ------------------------------------------------
            # Monta DataFrame
            # ------------------------------------------------

            linhas = []

            for r in resultados:

                if r["xs"] is None:

                    linhas.append({

                        "Pb (kPa)": r["Pb"] / 1000,

                        "xs (m)": np.nan,

                        "M1": np.nan,

                        "M2": np.nan,

                        "P1 (kPa)": np.nan,

                        "P2 (kPa)": np.nan,

                        "Me": np.nan,

                        "Pe (kPa)": np.nan,

                        "Erro (Pa)": np.nan,

                        "Erro (%)": np.nan

                    })

                    continue


                linhas.append({

                    "Pb (kPa)": r["Pb"] / 1000,

                    "xs (m)": r["xs"],

                    "M1": r["M1"],

                    "M2": r["M2"],

                    "P1 (kPa)": r["P1"] / 1000,

                    "P2 (kPa)": r["P2"] / 1000,

                    "Me": r["Me"],

                    "Pe (kPa)": r["Pe"] / 1000,

                    "Erro (Pa)": r["erro_abs"],

                    "Erro (%)": r["erro_rel"]

                })


            df = pd.DataFrame(linhas)


            st.subheader(
                "Tabela de validação"
            )

            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True
            )


            # ------------------------------------------------
            # Download CSV
            # ------------------------------------------------

            csv = df.to_csv(
                index=False
            ).encode("utf-8")


            st.download_button(
                label="EXPORTAR CSV",
                data=csv,
                file_name="validacao_bocal.csv",
                mime="text/csv",
                use_container_width=True
            )


            # ------------------------------------------------
            # Dados válidos
            # ------------------------------------------------

            df_validos = df.dropna(
                subset=["xs (m)"]
            )


            if len(df_validos) > 0:


                # --------------------------------------------
                # Gráfico Pb x xs
                # --------------------------------------------

                st.subheader(
                    "Posição do choque × Contrapressão"
                )


                fig_xs, ax_xs = plt.subplots(
                    figsize=(10, 5)
                )


                ax_xs.plot(
                    df_validos["Pb (kPa)"],
                    df_validos["xs (m)"],
                    marker="o",
                    linewidth=2
                )


                ax_xs.set_xlabel(
                    "Contrapressão Pb (kPa)"
                )

                ax_xs.set_ylabel(
                    "Posição do choque xs (m)"
                )

                ax_xs.set_title(
                    "Posição da onda de choque em função da contrapressão"
                )

                ax_xs.grid(
                    alpha=0.3
                )


                st.pyplot(
                    fig_xs,
                    use_container_width=True
                )


                plt.close(fig_xs)


                # --------------------------------------------
                # Gráfico do erro
                # --------------------------------------------

                st.subheader(
                    "Erro relativo"
                )


                fig_erro, ax_erro = plt.subplots(
                    figsize=(10, 5)
                )


                ax_erro.plot(
                    df_validos["Pb (kPa)"],
                    df_validos["Erro (%)"],
                    marker="o",
                    linewidth=2
                )


                ax_erro.set_xlabel(
                    "Contrapressão Pb (kPa)"
                )

                ax_erro.set_ylabel(
                    "Erro relativo (%)"
                )

                ax_erro.set_title(
                    "Erro relativo da condição Pe = Pb"
                )

                ax_erro.grid(
                    alpha=0.3
                )


                st.pyplot(
                    fig_erro,
                    use_container_width=True
                )


                plt.close(fig_erro)


                # --------------------------------------------
                # Estatísticas
                # --------------------------------------------

                erro_max = df_validos[
                    "Erro (%)"
                ].abs().max()


                erro_medio = df_validos[
                    "Erro (%)"
                ].abs().mean()


                col1, col2 = st.columns(2)


                with col1:

                    st.metric(
                        "Maior erro relativo",
                        f"{erro_max:.8f} %"
                    )


                with col2:

                    st.metric(
                        "Erro relativo médio",
                        f"{erro_medio:.8f} %"
                    )