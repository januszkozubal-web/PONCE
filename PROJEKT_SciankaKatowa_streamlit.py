import streamlit as st
import numpy as np

from PROJEKT_SciankaKatowa import oblicz_i_rysuj


def main():
    st.set_page_config(
        page_title="Ścianka oporowa – parcie czynne",
        layout="wide",
    )

    st.title("Ścianka oporowa – parcie czynne (Rankine, Poncelet/Coulomb)")

    st.markdown(
        """
        Wybierz, czy chcesz policzyć **parcie Rankine'a** czy **Ponceleta (Coulomba)**.
        Dzięki temu wykres jest większy i czytelniejszy.
        """
    )

    tryb = st.radio(
        "Tryb obliczeń",
        options=["Rankine", "Poncelet (Coulomb)"],
        index=0,
        horizontal=True,
    )

    with st.form("dane_wejsciowe"):
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Dane wspólne")
            q = st.number_input("q [kN/m²]", value=10.0, step=0.5)
            gamma = st.number_input("γ [kN/m³]", value=18.0, step=0.5)
            phi = st.slider("φ [°]", min_value=10.0, max_value=50.0, value=30.0, step=1.0)
            epsilon = st.slider("ε [°] (nachylenie naziomu)", min_value=-5.0, max_value=15.0, value=0.0, step=0.5)

            gamma_Q = st.selectbox(
                "γ_Q (wsp. bezp. obciążenia)",
                options=[1.0, 1.11, 1.2, 1.35, 1.5],
                index=0,
            )
            gamma_G = st.selectbox(
                "γ_G (wsp. bezp. gruntu)",
                options=[1.0, 1.2, 1.35, 1.5],
                index=0,
            )

            z_max = st.number_input("z_max [m] (wysokość ściany dla Rankine'a)", value=5.0, step=0.5, min_value=0.5)

        with col2:
            st.subheader("Poncelet (Coulomb) – potrzebne tylko, gdy wybrano Poncelet")
            beta = st.slider("β [°] (nachylenie ściany od pionu)", min_value=-30.0, max_value=30.0, value=10.0, step=0.5)

            delta_label = st.selectbox(
                "δ (tarcie grunt–ściana)",
                options=["0°", "φ/3", "2φ/3", "φ"],
                index=2,
            )

            l_max = st.number_input("l_max [m] (wysokość ściany dla Ponceleta)", value=5.0, step=0.5, min_value=0.5)

        submitted = st.form_submit_button("Oblicz")

    if submitted:
        DELTA_FRAC_MAP = {"0°": 0, "φ/3": 1, "2φ/3": 2, "φ": 3}

        try:
            d = {
                "q": float(q),
                "gamma": float(gamma),
                "phi_st": float(phi),
                "epsilon_st": float(epsilon),
                "gamma_Q": float(gamma_Q),
                "gamma_G": float(gamma_G),
                "z_max": float(z_max),
                "beta_st": float(beta),
                "delta_frac": DELTA_FRAC_MAP[delta_label],
                "l_max": float(l_max),
            }

            path1, path2, wyniki, fig1, fig2 = oblicz_i_rysuj(d, return_figures=True)

            st.success("Obliczenia zakończone pomyślnie.")

            # powiększenie wykresów na potrzeby wyświetlenia w Streamlit
            if tryb == "Rankine":
                fig1.set_size_inches(10.0, 14.0)
                # w przeglądarce większa wysokość – bez wymuszania proporcji 1:1
                if fig1.axes:
                    fig1.axes[0].set_aspect("auto")
            else:  # "Poncelet (Coulomb)"
                fig2.set_size_inches(10.0, 14.0)
                if fig2.axes:
                    fig2.axes[0].set_aspect("auto")

            if tryb == "Rankine":
                st.subheader("Wyniki – Rankine")
                r = wyniki["rankine"]
                st.write(f"E_a,d = {r['E_a']:.2f} kN/m")
                st.write(f"E_ah,d = {r['E_ah']:.2f} kN/m")
                st.write(f"E_av,d = {r['E_av']:.2f} kN/m")
                st.write(f"z_c = {r['z_c']:.2f} m (od góry)")
                if "K_a" in r:
                    st.write(f"K_a = {r['K_a']:.3f}")
                if "omega_epsilon_deg" in r:
                    st.write(f"ω_ε = {r['omega_epsilon_deg']:.2f}°")
                    st.write(f"θ = {r['theta_deg']:.2f}°")
                    st.write(f"90° − θ = {r['comp_90_minus_theta_deg']:.2f}°")
                st.pyplot(fig1, clear_figure=False)

            else:  # "Poncelet (Coulomb)"
                st.subheader("Wyniki – Poncelet (Coulomb)")
                p = wyniki["poncelet"]
                st.write(f"E_a,d = {p['E_a']:.2f} kN/m")
                st.write(f"E_ah,d = {p['E_ah']:.2f} kN/m")
                st.write(f"E_av,d = {p['E_av']:.2f} kN/m")
                st.write(f"z_c = {p['z_c']:.2f} m (od góry)")
                # Dla Ponceleta przydatne są także kąty z części Rankine'a
                r = wyniki.get("rankine", {})
                if "omega_epsilon_deg" in r:
                    st.write(f"ω_ε = {r['omega_epsilon_deg']:.2f}°")
                    st.write(f"θ = {r['theta_deg']:.2f}°")
                    st.write(f"90° − θ = {r['comp_90_minus_theta_deg']:.2f}°")
                st.pyplot(fig2, clear_figure=False)

            st.caption(f"Pliki PDF zostały zapisane jako: {path1} oraz {path2}")

        except Exception as e:
            st.error(f"Błąd obliczeń: {e}")


if __name__ == "__main__":
    main()

