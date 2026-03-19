# -*- coding: utf-8 -*-
# =============================================================================
# PROJEKT_SciankaKatowa.py
# Scianka oporowa - konstrukcje betonowe, temat 5
# Parcie czynne: Rankine oraz Coulomb / Poncelet
# Uruchom: python PROJEKT_SciankaKatowa.py  -> od razu okienko, potem dwie grafiki PDF
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KATALOG = os.path.dirname(os.path.abspath(__file__))


def oblicz_i_rysuj(d, return_figures=False, save_pdf=True):
    """Obliczenia Rankine + Poncelet i (opcjonalnie) zapis dwóch PDF.

    Jeśli `return_figures=True`, zwraca też `fig1`, `fig2` (bez zamykania).
    Jeśli `save_pdf=False`, nie wykonuje się zapisu do plików, a `path1/path2` ustawiane są na `None`.
    """
    pi = np.pi
    q = d["q"]
    gamma = d["gamma"]
    phi = d["phi_st"] * pi / 180
    epsilon = d["epsilon_st"] * pi / 180
    gamma_Q = d["gamma_Q"]
    gamma_G = d["gamma_G"]
    z_max = d["z_max"]
    beta = d["beta_st"] * pi / 180
    delta_frac = d["delta_frac"]
    l_max = d["l_max"]

    q_d = q * gamma_Q
    gamma_d = gamma * gamma_G

    # Wspólny współczynnik do zmniejszenia długości strzałek / wektorów (wizualnie ~5× mniejsze)
    VIS_SCALE = 0.2
    # Osobny współczynnik dla obwiedni naprężeń e_a (wizualnie ~5× mniejsze)
    EA_VIS_SCALE = 0.2

    # ----- Rankine -----
    # Uogólniona postać z kątem ε; dla ε → 0 stosujemy klasyczny wzór Rankine'a
    # tylko do obliczenia K_a_gamma oraz h_z_d (bez kłopotów numerycznych),
    # natomiast kąt ω_ε liczymy zawsze z ogólnego wzoru.
    omega_epsilon = np.arcsin(np.sin(epsilon) / np.sin(phi))

    if abs(epsilon) < 1e-6:
        # klasyczne parcie czynne Rankine'a przy poziomym naziomie
        K_a_gamma = np.tan(np.pi / 4 - phi / 2) ** 2
        h_z_d = q_d / gamma_d
    else:
        K_a_gamma = np.sin(omega_epsilon - epsilon) * np.cos(epsilon) / np.sin(omega_epsilon + epsilon)
        h_z_d = q_d / (gamma_d * np.cos(epsilon))

    # Kąt θ potrzebny przy konstrukcji klinu (Poncelet):
    #   θ = π/4 + φ/2 + (ω_ε − ε)/2
    theta = pi / 4 + phi / 2 + (omega_epsilon - epsilon) / 2

    # Wypisanie ω_ε, θ oraz 90° − θ (w stopniach)
    omega_deg = np.degrees(omega_epsilon)
    theta_deg = np.degrees(theta)
    print(f"omega_epsilon = {omega_deg:.2f}°")
    print(f"theta = {theta_deg:.2f}°")
    print(f"90° - theta = {90.0 - theta_deg:.2f}°")

    # Wypisanie współczynnika parcia czynnego Rankine'a
    print(f"K_a (Rankine, od ciężaru gruntu) = {K_a_gamma:.4f}")

    def e_a(z):
        return K_a_gamma * (z + h_z_d) * gamma_d

    E_a = K_a_gamma * gamma_d * (z_max**2 / 2 + h_z_d * z_max)
    E_ah = E_a * np.cos(epsilon)
    E_av = E_a * np.sin(epsilon)  # składowa pionowa
    moment_1 = K_a_gamma * gamma_d * (z_max**3 / 3 + h_z_d * z_max**2 / 2)
    z_c = moment_1 / E_a

    # ----- Poncelet (Coulomb) -----
    delta_opt = np.array([0, 1/3, 2/3, 1])
    delta = phi * delta_opt[delta_frac]

    D_coulomb = np.sin(phi + delta) * np.sin(phi - epsilon) / (
        np.cos(beta + delta) * np.cos(epsilon - beta)
    )
    if D_coulomb < 0:
        raise ValueError("Coulomb: D < 0 (sprawdź kąty, np. epsilon < phi)")

    K_agamma_coulomb = np.cos(phi - beta)**2 / (
        np.cos(beta + delta) * (1 + np.sqrt(D_coulomb))**2
    )
    K_aq_coulomb = K_agamma_coulomb / np.cos(epsilon - beta)

    # Wypisanie współczynników parcia czynnego Coulomba / Ponceleta
    print(f"K_a,γ (Coulomb) = {K_agamma_coulomb:.4f}")
    print(f"K_a,q (Coulomb) = {K_aq_coulomb:.4f}")

    def e_a_poncelet(l):
        return q_d * K_aq_coulomb + K_agamma_coulomb * gamma_d * l

    E_a_poncelet = q_d * K_aq_coulomb * l_max + K_agamma_coulomb * gamma_d * (l_max**2 / 2)
    moment_poncelet = q_d * K_aq_coulomb * (l_max**2 / 2) + K_agamma_coulomb * gamma_d * (l_max**3 / 3)
    z_c_poncelet = moment_poncelet / E_a_poncelet
    E_ah_poncelet = E_a_poncelet * np.cos(beta + delta)
    E_av_poncelet = E_a_poncelet * np.sin(beta + delta)  # składowa pionowa

    # ----- Wykres Rankine -----
    z_plot = np.linspace(0, z_max, 200)
    ea_plot = e_a(z_plot)
    ea_plot_vis = EA_VIS_SCALE * ea_plot
    x_wall_r = np.zeros_like(z_plot)
    y_wall_r = z_plot
    # Do rysowania używamy przeskalowanego naprężenia (nie wpływa to na E_a itd.)
    x_ea_r = x_wall_r + ea_plot_vis * np.cos(epsilon)
    y_ea_r = y_wall_r - ea_plot_vis * np.sin(epsilon)
    # Odporne na NaN/Inf wyznaczanie x_max
    finite_x_ea_r = x_ea_r[np.isfinite(x_ea_r)]
    if finite_x_ea_r.size == 0:
        x_max = 1.0
    else:
        x_max = np.max(finite_x_ea_r) * 1.25
    # Marginesy: duży zapas od góry (żeby rysunek nie był ucięty)
    y_hi = -3.0
    y_lo = z_max + 1.0

    fig1, ax1 = plt.subplots(figsize=(9, 7))
    ax1.set_position([0.08, 0.04, 0.88, 0.90])  # osie niżej – minimalny odstęp od założeń
    # mniejszy margines z lewej strony (bliżej do końca strzałki E_ah)
    ax1.set_xlim(-x_max * 0.4, x_max)
    ax1.set_ylim(z_max + 1.0, y_hi)  # grunt 1 m poniżej podstawy ściany
    ax1.set_xlabel("")  # brak opisu osi poziomej
    ax1.set_ylabel("głębokość z [m]")
    ax1.set_title("Parcie czynne Rankine'a")
    ax1.set_aspect("equal")

    # Rysunek naziomu (ε) nad ścianą – po stronie parcia
    ground_len = x_max * 0.5
    x_ground = np.array([0.0, ground_len])
    y_ground = np.array([0.0, -np.tan(epsilon) * ground_len])
    ax1.plot(x_ground, y_ground, color="saddlebrown", lw=2)

    # Obszar gruntu – żółty, między naziomem a poziomem 1 m poniżej podstawy ściany
    y_soil_bottom = z_max + 1.0
    x_soil = np.array([0.0, ground_len, ground_len, 0.0])
    y_soil = np.array([0.0, y_ground[1], y_soil_bottom, y_soil_bottom])
    ax1.fill(x_soil, y_soil, color="#ffdd77", alpha=0.7, edgecolor="none")

    # Parcie – szare
    ax1.fill(
        np.concatenate([x_wall_r, x_ea_r[::-1]]),
        np.concatenate([y_wall_r, y_ea_r[::-1]]),
        color="#bbbbbb",
        alpha=0.7,
        edgecolor="none",
    )
    ax1.plot(x_wall_r, y_wall_r, lw=3, color="#333333", label="Ściana")
    ax1.plot(x_ea_r, y_ea_r, lw=2.5, color="#666666")
    ax1.plot([x_ea_r[0], x_wall_r[0]], [y_ea_r[0], y_wall_r[0]], color="#666666", lw=1)
    ax1.plot([x_ea_r[-1], x_wall_r[-1]], [y_ea_r[-1], y_wall_r[-1]], color="#666666", lw=1)

    # Strzałki q – obciążenie naziomu, gęstsze, długość proporcjonalna do q (w dół)
    n_arrows = 8
    xs_ar = np.linspace(0.05 * ground_len, 0.95 * ground_len, n_arrows)
    ys_ar = -np.tan(epsilon) * xs_ar
    # Skala długości strzałek na podstawie q_d i sił z wykresu
    scale_q = VIS_SCALE * 0.25 * x_max * (q_d / max(E_ah, 1e-6))
    if scale_q > 1e-8:
        for xa, ya in zip(xs_ar, ys_ar):
            # strzałka od punktu powyżej naziomu w dół do naziomu (grot na naziomie)
            ax1.annotate(
                "",
                xy=(xa, ya),                # grot na naziomie
                xytext=(xa, ya - scale_q),  # początek powyżej naziomu
                arrowprops=dict(arrowstyle="->", color="black", lw=1.1),
            )

    # Oznaczenie kąta ε (odchylenie naziomu od poziomu)
    epsilon_deg = np.degrees(epsilon)
    arc_r = 0.9
    theta1 = 0.0
    theta2 = epsilon_deg
    arc_eps = matplotlib.patches.Arc(
        (0.0, 0.0),
        width=2 * arc_r,
        height=2 * arc_r,
        angle=0.0,
        theta1=min(theta1, theta2),
        theta2=max(theta1, theta2),
        color="black",
        lw=1.0,
    )
    ax1.add_patch(arc_eps)
    ax1.text(
        0.6 * arc_r,
        0.25 * arc_r,
        "ε",
        fontsize=11,
        ha="center",
        va="center",
    )

    # Opis wartości parcia w górze i u dołu ściany (Rankine)
    ea_top = ea_plot[0]
    ea_bot = ea_plot[-1]
    # współrzędne na wykresie parcia (na obwiedni)
    x_top = x_ea_r[0]
    y_top = y_ea_r[0]
    x_bot = x_ea_r[-1]
    y_bot = y_ea_r[-1]
    # lekkie odsunięcie tekstów od obwiedni
    dx = 0.03 * x_max
    ax1.text(
        x_top + dx,
        y_top - 0.2,
        f"e_a, górą = {ea_top:.2f} kN/m²",
        fontsize=8,
        ha="left",
        va="top",
        color="#333333",
    )
    ax1.text(
        x_bot + dx,
        y_bot + 0.2,
        f"e_a, dołem = {ea_bot:.2f} kN/m²",
        fontsize=8,
        ha="left",
        va="bottom",
        color="#333333",
    )

    # Skala siły: ta sama dla poziomej i pionowej (np. 0.25*x_max na długość E_ah)
    scale_E = VIS_SCALE * 0.25 * x_max / E_ah
    arrow_x_end = -E_ah * scale_E
    ax1.annotate("", xy=(arrow_x_end, z_c), xytext=(0, z_c),
                 arrowprops=dict(arrowstyle="->", color="darkred", lw=2.5))
    # E_ah,d w górę nad strzałką (nie w bok)
    ax1.text(0.5 * arrow_x_end, z_c + 0.95, f"$E_{{ah,d}}$ = {E_ah:.2f} kN/m", fontsize=9, color="darkred", fontweight="bold", va="bottom", ha="center")
    # Strzałka E_av,d dokładnie w dół, w tej samej skali co pozioma
    scale_av = E_av * scale_E
    z_av_tip = z_c + scale_av
    ax1.annotate("", xy=(0, z_av_tip), xytext=(0, z_c),
                 arrowprops=dict(arrowstyle="->", color="darkred", lw=2.5))
    ax1.text(0.12, z_c + scale_av / 2, f"$E_{{av,d}}$ = {E_av:.2f} kN/m", fontsize=9, color="darkred", va="center")
    ax1.plot(0, z_c, "o", color="darkred", markersize=8)
    # Blok wartości w prawym dolnym rogu (z dopisanym współczynnikiem K_a)
    txt_rankine = (
        f"$E_{{a,d}}$ = {E_a:.2f} kN/m\n"
        f"$E_{{ah,d}}$ = {E_ah:.2f} kN/m\n"
        f"$E_{{av,d}}$ = {E_av:.2f} kN/m\n"
        f"$z_c$ = {z_c:.2f} m (od góry)\n"
        f"$h_{{z,d}}$ = {h_z_d:.2f} m\n"
        f"$K_a$ = {K_a_gamma:.3f}"
    )
    ax1.text(
        0.98,
        0.02,
        txt_rankine,
        transform=ax1.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="right",
        color="#222222",
    )
    path1 = os.path.join(KATALOG, "PROJEKT_SciankaKatowa_wykres_parcia.pdf")
    if save_pdf:
        fig1.savefig(path1, format="pdf", bbox_inches="tight")
    else:
        path1 = None
    if not return_figures:
        plt.close(fig1)

    # ----- Wykres Poncelet -----
    l_plot = np.linspace(0, l_max, 200)
    ea_p_plot = e_a_poncelet(l_plot)
    ea_p_plot_vis = EA_VIS_SCALE * ea_p_plot
    x_wall_p = l_plot * np.sin(beta)
    y_wall_p = l_plot * np.cos(beta)
    # Do rysowania używamy przeskalowanej obwiedni naprężeń
    x_ea_p = x_wall_p + ea_p_plot_vis * np.cos(beta + delta)
    y_ea_p = y_wall_p - ea_p_plot_vis * np.sin(beta + delta)
    finite_x_ea_p = x_ea_p[np.isfinite(x_ea_p)]
    if finite_x_ea_p.size == 0:
        x_max_p = 1.0
    else:
        x_max_p = np.max(finite_x_ea_p) * 1.25
    # Zakres y z zawartości rysunku – cały wykres widoczny (góra nie ucięta)
    y_min_content = min(0.0, np.min(y_wall_p), np.min(y_ea_p))
    y_max_content = max(np.max(y_wall_p), np.max(y_ea_p))
    margin = 1.2
    y_hi_p = y_min_content - margin
    y_lo_p = y_max_content + margin

    fig2, ax2 = plt.subplots(figsize=(9, 7))
    ax2.set_position([0.08, 0.04, 0.88, 0.90])
    # mniejszy margines z lewej strony
    ax2.set_xlim(-x_max_p * 0.4, x_max_p)
    ax2.set_ylim(y_lo_p, y_hi_p)
    ax2.set_xlabel("")  # brak opisu osi poziomej
    ax2.set_ylabel("głębokość rzutu l·cos(β) [m]")
    ax2.set_title("Parcie czynne Ponceleta (Coulomb)")
    ax2.set_aspect("equal")

    # Naziom (ε) dla Ponceleta – po stronie parcia
    ground_len_p = x_max_p * 0.5
    x_ground_p = np.array([0.0, ground_len_p])
    y_ground_p = np.array([0.0, -np.tan(epsilon) * ground_len_p])
    ax2.plot(x_ground_p, y_ground_p, color="saddlebrown", lw=2)

    # Obszar gruntu – żółty, ograniczony:
    #  - ścianą,
    #  - naziomem,
    #  - pionową linią pod końcem naziomu,
    #  - poziomą linią u dołu (ok. 1 m poniżej stopy ściany),
    #  - oraz fragmentem poziomym pod ścianą.
    x_wall_foot = x_wall_p[-1]
    y_wall_foot = y_wall_p[-1]
    y_soil_bottom_p = y_wall_foot + 1.0  # ~1 m poniżej stopy ściany
    x_soil_p = np.array([
        0.0,                 # wierzchołek ściany / początek naziomu
        x_wall_foot,         # wzdłuż ściany do stopy
        x_wall_foot,         # pionowo w dół
        x_ground_p[1],       # poziomo pod naziomem
        x_ground_p[1],       # pionowo w górę do końca naziomu
    ])
    y_soil_p = np.array([
        0.0,                 # wierzchołek ściany
        y_wall_foot,         # stopa ściany
        y_soil_bottom_p,     # dół przy stopie ściany
        y_soil_bottom_p,     # dół przy pionie pod końcem naziomu
        y_ground_p[1],       # koniec naziomu
    ])
    ax2.fill(x_soil_p, y_soil_p, color="#ffdd77", alpha=0.7, edgecolor="none")

    # Parcie – zachowaj zielony klin
    ax2.fill(
        np.concatenate([x_wall_p, x_ea_p[::-1]]),
        np.concatenate([y_wall_p, y_ea_p[::-1]]),
        color="darkgreen",
        alpha=0.35,
        edgecolor="none",
    )
    ax2.plot(x_wall_p, y_wall_p, lw=3, color="#333333")
    ax2.plot(x_ea_p, y_ea_p, lw=2.5, color="darkgreen")
    ax2.plot([x_ea_p[0], x_wall_p[0]], [y_ea_p[0], y_wall_p[0]], color="darkgreen", lw=1)
    ax2.plot([x_ea_p[-1], x_wall_p[-1]], [y_ea_p[-1], y_wall_p[-1]], color="darkgreen", lw=1)
    x_c_wall = z_c_poncelet * np.sin(beta)
    y_c_wall = z_c_poncelet * np.cos(beta)

    # Strzałki q – obciążenie naziomu (gęściej, długość proporcjonalna do q, w dół)
    n_arrows_p = 8
    xs_ar_p = np.linspace(0.05 * ground_len_p, 0.95 * ground_len_p, n_arrows_p)
    ys_ar_p = -np.tan(epsilon) * xs_ar_p
    scale_q_p = VIS_SCALE * 0.25 * x_max_p * (q_d / max(E_ah_poncelet, 1e-6))
    if scale_q_p > 1e-8:
        for xa, ya in zip(xs_ar_p, ys_ar_p):
            # strzałka od punktu powyżej naziomu w dół do naziomu (grot na naziomie)
            ax2.annotate(
                "",
                xy=(xa, ya),                   # grot na naziomie
                xytext=(xa, ya - scale_q_p),   # początek powyżej naziomu
                arrowprops=dict(arrowstyle="->", color="black", lw=1.1),
            )

    # Dodatkowa linia równoległa do naziomu, zaczynająca się od dołu ściany,
    # obciążona stałą wartością q + γ * l_max * cos(β)
    x_wall_foot = x_wall_p[-1]
    y_wall_foot = y_wall_p[-1]
    x_ground2_p = np.array([x_wall_foot, x_wall_foot + ground_len_p])
    y_ground2_p = np.array([y_wall_foot, y_wall_foot - np.tan(epsilon) * ground_len_p])
    ax2.plot(x_ground2_p, y_ground2_p, color="saddlebrown", lw=1.5, linestyle="--")

    # Strzałki na tej linii – stałe: q_eff = q_d + γ_d * l_max * cos(β)
    n_arrows2_p = 8
    xs2 = np.linspace(x_ground2_p[0] + 0.05 * ground_len_p,
                      x_ground2_p[0] + 0.95 * ground_len_p,
                      n_arrows2_p)
    ys2 = y_wall_foot - np.tan(epsilon) * (xs2 - x_wall_foot)
    q_eff = q_d + gamma_d * l_max * np.cos(beta)
    # długość strzałki proporcjonalna do q_eff (skalowana względem pierwotnego q_d)
    scale_q_eff = scale_q_p * (q_eff / max(q_d, 1e-6)) if q_d != 0 else 0.0
    if scale_q_eff > 1e-8:
        for xa, ya in zip(xs2, ys2):
            ax2.annotate(
                "",
                xy=(xa, ya),                        # grot na linii równoległej do naziomu
                xytext=(xa, ya - scale_q_eff),      # początek powyżej linii
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
            )

    # Wartość efektywnego obciążenia q_eff na dodatkowej linii
    ax2.text(
        x_ground2_p[1] + 0.05 * x_max_p,
        y_ground2_p[1],
        f"q' = {q_eff:.2f} kN/m²",
        fontsize=8,
        ha="left",
        va="center",
        color="#222222",
    )

    # Opis wartości parcia w górze i u dołu ściany (Poncelet)
    ea_p_top = ea_p_plot[0]
    ea_p_bot = ea_p_plot[-1]
    x_p_top = x_ea_p[0]
    y_p_top = y_ea_p[0]
    x_p_bot = x_ea_p[-1]
    y_p_bot = y_ea_p[-1]
    dx_p = 0.03 * x_max_p
    ax2.text(
        x_p_top + dx_p,
        y_p_top - 0.2,
        f"e_a, górą = {ea_p_top:.2f} kN/m²",
        fontsize=8,
        ha="left",
        va="top",
        color="#222222",
    )
    ax2.text(
        x_p_bot + dx_p,
        y_p_bot + 0.2,
        f"e_a, dołem = {ea_p_bot:.2f} kN/m²",
        fontsize=8,
        ha="left",
        va="bottom",
        color="#222222",
    )

    # Kąt β – odchylenie ściany od pionu (przy wierzchołku ściany)
    beta_deg = np.degrees(beta)
    arc_r_b = 0.9
    theta_v = 90.0
    theta_wall = 90.0 - beta_deg
    theta1_b = min(theta_v, theta_wall)
    theta2_b = max(theta_v, theta_wall)
    arc_beta = matplotlib.patches.Arc(
        (0.0, 0.0),
        width=2 * arc_r_b,
        height=2 * arc_r_b,
        angle=0.0,
        theta1=theta1_b,
        theta2=theta2_b,
        color="black",
        lw=1.0,
    )
    ax2.add_patch(arc_beta)
    ax2.text(
        arc_r_b * np.sin(np.radians(beta_deg / 2.0)),
        arc_r_b * np.cos(np.radians(beta_deg / 2.0)),
        r"$\beta$",
        fontsize=11,
        ha="center",
        va="center",
    )

    # Oś wzdłuż ściany (Poncelet) – z opisem l = 0 i l = l_max
    # Oś kończy się ok. 0.5 m poniżej stopy ściany (nie za daleko)
    l_axis = z_max + 0.5
    x_axis_wall_top = 0.0
    y_axis_wall_top = 0.0
    x_axis_wall_bot = l_axis * np.sin(beta)
    y_axis_wall_bot = l_axis * np.cos(beta)
    # Strzałka skierowana ku górze (przeciwny zwrot niż wcześniej)
    ax2.annotate(
        "",
        xy=(x_axis_wall_bot, y_axis_wall_bot),
        xytext=(x_axis_wall_top, y_axis_wall_top),
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    ax2.text(
        x_axis_wall_top - 0.1,
        y_axis_wall_top - 0.1,
        "l = 0",
        fontsize=8,
        ha="right",
        va="top",
        color="#222222",
    )
    ax2.text(
        x_axis_wall_bot - 0.1,
        y_axis_wall_bot + 0.1,
        f"l = {l_max:.2f} m",
        fontsize=8,
        ha="right",
        va="bottom",
        color="#222222",
    )

    # Normalna do ściany i kąt δ (Poncelet) – zaznaczone przy środku działania E_a
    # Kierunek normalnej wychodzącej z gruntu (linia robocza – bez grotu)
    n_vec = np.array([-np.cos(beta), np.sin(beta)])
    n_len = 1.4
    x_n_end = x_c_wall + n_len * n_vec[0]
    y_n_end = y_c_wall + n_len * n_vec[1]
    ax2.plot([x_c_wall, x_n_end], [y_c_wall, y_n_end], color="black", lw=1.4)

    # Linia odchylenia normalnej o kąt δ – cienka linia w kierunku parcia
    d_len = 1.1
    x_d_end = x_c_wall + d_len * np.cos(beta + delta)
    y_d_end = y_c_wall - d_len * np.sin(beta + delta)
    ax2.plot([x_c_wall, x_d_end], [y_c_wall, y_d_end], color="black", lw=1.0, linestyle="--")

    # Kierunek parcia (wzdłuż E_a – β+δ od poziomu) – strzałka na linii odchylenia
    r_len = 0.8
    x_r_end = x_c_wall + r_len * np.cos(beta + delta)
    y_r_end = y_c_wall - r_len * np.sin(beta + delta)
    ax2.annotate(
        "",
        xy=(x_r_end, y_r_end),
        xytext=(x_c_wall, y_c_wall),
        arrowprops=dict(arrowstyle="->", color="darkred", lw=1.6),
    )

    # Łuk dla kąta δ pomiędzy normalną a kierunkiem parcia
    ang_n = np.degrees(np.arctan2(n_vec[1], n_vec[0]))
    ang_r = np.degrees(np.arctan2(-np.sin(beta + delta), np.cos(beta + delta)))
    theta1_d = min(ang_n, ang_r)
    theta2_d = max(ang_n, ang_r)
    arc_r_d = 0.6
    arc_delta = matplotlib.patches.Arc(
        (x_c_wall, y_c_wall),
        width=2 * arc_r_d,
        height=2 * arc_r_d,
        angle=0.0,
        theta1=theta1_d,
        theta2=theta2_d,
        color="black",
        lw=1.0,
    )
    ax2.add_patch(arc_delta)
    # Pozycja etykiety δ – w środku łuku
    ang_mid = 0.5 * (theta1_d + theta2_d)
    x_delta_txt = x_c_wall + arc_r_d * np.cos(np.radians(ang_mid))
    y_delta_txt = y_c_wall + arc_r_d * np.sin(np.radians(ang_mid))
    ax2.text(
        x_delta_txt,
        y_delta_txt,
        r"$\delta$",
        fontsize=11,
        ha="center",
        va="center",
    )

    # Kąt ε – odchylenie naziomu od poziomu, zaznaczony osobno przy końcu naziomu
    epsilon_deg_p = np.degrees(epsilon)
    arc_r_eps_p = 0.7
    # łuk między poziomem a linią naziomu (po prawej stronie rysunku)
    theta1_eps_p = 0.0
    theta2_eps_p = -epsilon_deg_p
    arc_eps_p = matplotlib.patches.Arc(
        (x_ground_p[1], y_ground_p[1]),
        width=2 * arc_r_eps_p,
        height=2 * arc_r_eps_p,
        angle=0.0,
        theta1=min(theta1_eps_p, theta2_eps_p),
        theta2=max(theta1_eps_p, theta2_eps_p),
        color="black",
        lw=1.0,
    )
    ax2.add_patch(arc_eps_p)
    ax2.text(
        x_ground_p[1] - 0.6 * arc_r_eps_p,
        y_ground_p[1] - 0.25 * arc_r_eps_p,
        "ε",
        fontsize=11,
        ha="center",
        va="center",
    )

    # Ta sama skala dla E_ah i E_av (jeśli trzeba skrócić, skracamy OBA wektory, żeby zachować proporcje)
    scale_E_p = VIS_SCALE * 0.25 * x_max_p / E_ah_poncelet
    full_len_h = E_ah_poncelet * scale_E_p
    full_len_v = E_av_poncelet * scale_E_p

    # maksymalna dopuszczalna długość w pionie (żeby strzałka nie wyszła poza rysunek)
    max_len_v = max(0, (y_lo_p - y_c_wall) * 0.85)
    if max_len_v > 0 and abs(full_len_v) > max_len_v:
        factor = max_len_v / abs(full_len_v)
        full_len_v *= factor
        full_len_h *= factor

    arrow_x_end_p = x_c_wall - full_len_h
    ax2.annotate("", xy=(arrow_x_end_p, y_c_wall), xytext=(x_c_wall, y_c_wall),
                 arrowprops=dict(arrowstyle="->", color="darkred", lw=2.5))
    # E_ah,d w górę nad strzałką (nie w bok)
    x_mid_p = 0.5 * (x_c_wall + arrow_x_end_p)
    ax2.text(x_mid_p, y_c_wall + 0.95, f"$E_{{ah,d}}$ = {E_ah_poncelet:.2f} kN/m", fontsize=9, color="darkred", fontweight="bold", va="bottom", ha="center")
    # Strzałka E_av,d w dół (wzdłuż ściany w głąb); ta sama skala co dla poziomej
    delta_deg = np.degrees(delta)
    y_av_tip = y_c_wall + full_len_v
    ax2.annotate("", xy=(x_c_wall, y_av_tip), xytext=(x_c_wall, y_c_wall),
                 arrowprops=dict(arrowstyle="->", color="darkred", lw=2.5))
    ax2.text(x_c_wall + 0.15, (y_c_wall + y_av_tip) / 2, f"$E_{{av,d}}$ = {E_av_poncelet:.2f} kN/m", fontsize=9, color="darkred", va="center")
    ax2.plot(x_c_wall, y_c_wall, "o", color="darkred", markersize=8)
    # Blok wartości (z dopisanymi współczynnikami K_a)
    txt_poncelet = (
        f"$E_{{a,d}}$ = {E_a_poncelet:.2f} kN/m\n"
        f"$E_{{ah,d}}$ = {E_ah_poncelet:.2f} kN/m\n"
        f"$E_{{av,d}}$ = {E_av_poncelet:.2f} kN/m\n"
        f"$z_c$ = {z_c_poncelet:.2f} m (od góry)\n"
        f"$K_{{a,\\gamma}}$ = {K_agamma_coulomb:.3f}\n"
        f"$K_{{a,q}}$ = {K_aq_coulomb:.3f}"
    )
    ax2.text(
        0.98,
        0.02,
        txt_poncelet,
        transform=ax2.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        horizontalalignment="right",
        color="#222222",
    )
    path2 = os.path.join(KATALOG, "PROJEKT_SciankaKatowa_wykres_parcia_poncelet.pdf")
    if save_pdf:
        fig2.savefig(path2, format="pdf", bbox_inches="tight")
    else:
        path2 = None
    if not return_figures:
        plt.close(fig2)

    wyniki = {
        "rankine": {
            "E_a": E_a,
            "E_ah": E_ah,
            "E_av": E_av,
            "z_c": z_c,
            "h_z": h_z_d,
            "K_a": K_a_gamma,
            "omega_epsilon_deg": omega_deg,
            "theta_deg": theta_deg,
            "comp_90_minus_theta_deg": 90.0 - theta_deg,
        },
        "poncelet": {
            "E_a": E_a_poncelet,
            "E_ah": E_ah_poncelet,
            "E_av": E_av_poncelet,
            "z_c": z_c_poncelet,
        },
    }
    if return_figures:
        return path1, path2, wyniki, fig1, fig2
    return path1, path2, wyniki


def main():
    """Prosty tryb CLI: wywołaj obliczenia z domyślnym zestawem parametrów i zapisz PDF."""
    d = {
        "q": 10.0,
        "gamma": 18.0,
        "phi_st": 30.0,
        "epsilon_st": 0.0,
        "gamma_Q": 1.0,
        "gamma_G": 1.0,
        "z_max": 5.0,
        "beta_st": 10.0,
        "delta_frac": 2,  # 2φ/3
        "l_max": 5.0,
    }
    path1, path2, wyniki = oblicz_i_rysuj(d, return_figures=False, save_pdf=True)
    print("Zapisano pliki PDF:", path1, "oraz", path2)


if __name__ == "__main__":
    main()
