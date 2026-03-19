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
**Autor:** jvk &nbsp; | &nbsp; **Licencja:** MIT &nbsp; | &nbsp; **Wsparcie:** ćwiczeń projektowych

Funkcje napisane w Pythonie i `matplotlib` liczą parcie według założeń metody Rankine'a.
Tutaj ściana oporowa jest pionowa, a mamy do czynienia z nachylonym naziomem (kąt `ε`).
Powstaje wypadkowa nachylona pod tym samym kątem co naziom, dlatego wypadkowa ma składową poziomą i pionową.
        """
    )

    st.markdown(
        """
        Wybierz, czy chcesz policzyć **parcie Rankine'a** czy **Ponceleta (Coulomba)**.
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
            epsilon = st.slider("ε [°] (nachylenie naziomu)", min_value=-5.0, max_value=25.0, value=0.0, step=0.5)

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

    with st.expander("Wzory (Rankine / Poncelet)", expanded=False):
        # Sekcja opisowa - nie wykonuje obliczeń, tylko pokazuje wzory zgodne z `PROJEKT_SciankaKatowa.py`.
        if tryb == "Rankine":
            st.markdown(
                r"""
### Podstawowe wzory – Rankine

Niech:
`q_d = q * γ_Q`, `γ_d = γ * γ_G`, `φ` – kąt tarcia wewnętrznego, `ε` – nachylenie naziomu (od poziomu).

1. Kąt:
   `ω_ε = arcsin( sin(ε) / sin(φ) )`
2. Współczynnik parcia od ciężaru gruntu:
   `K_{a,γ} = sin(ω_ε − ε) cos(ε) / sin(ω_ε + ε)`  (dla `ε ≠ 0`)

   Dla `ε = 0` (naziom poziomy):
   `K_{a,γ} = tan^2(π/4 − φ/2)`
3. Wysokość charakterystyczna:
   `h_{z,d} = q_d / (γ_d cos(ε))` (dla `ε ≠ 0`), a dla `ε = 0`: `h_{z,d} = q_d / γ_d`
4. Jednostkowe parcie czynne w funkcji głębokości `z`:
   `e_a(z) = K_{a,γ} (z + h_{z,d}) γ_d`
                """
            )
        else:  # "Poncelet (Coulomb)"
            st.markdown(
                r"""
### Podstawowe wzory – Poncelet (Coulomb)

Metoda Ponceleta (klin poślizgu) prowadzi do postaci zbliżonej do Rankine'a, ale ze współczynnikami Coulomba.
Kąt `ω_ε` i `θ` liczymy z tej samej geometrii klinu.

Niech `β` oznacza odchylenie ściany od pionu (w kodzie `beta_st`), a `δ` – tarcie grunt–ściana.
`q_d = q * γ_Q`, `γ_d = γ * γ_G`.

1. Kąt:
   `ω_ε = arcsin( sin(ε) / sin(φ) )`
2. Kąt konstrukcji klina:
   `θ = π/4 + φ/2 + (ω_ε − ε)/2`
   `90° − θ = π/2 − θ`
3. Współczynniki parcia:
   `D = [sin(φ + δ) sin(φ − ε)] / [cos(β + δ) cos(ε − β)]`
   `K_{a,γ} = cos^2(φ − β) / ( cos(β + δ) (1 + sqrt(D))^2 )`
   `K_{a,q} = K_{a,γ} / cos(ε − β)`
4. Jednostkowe parcie czynne jako funkcja długości `l` wzdłuż ściany:
   `e_a(l) = q_d * K_{a,q} + K_{a,γ} * γ_d * l`
                """
            )

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

            path1, path2, wyniki, fig1, fig2 = oblicz_i_rysuj(d, return_figures=True, save_pdf=False)

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

            st.caption("Rysunki zostały wygenerowane (bez zapisu PDF).")

        except Exception as e:
            st.error(f"Błąd obliczeń: {e}")

if __name__ == "__main__":
    main()

