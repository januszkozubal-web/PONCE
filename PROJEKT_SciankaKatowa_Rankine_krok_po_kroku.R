# -----------------------------------------------------------------------------
# PROJEKT_SciankaKatowa_Rankine_krok_po_kroku.R
# Proste obliczenia Rankine'a krok po kroku (wariant jak w projekcie Python).
# Uruchom: Rscript PROJEKT_SciankaKatowa_Rankine_krok_po_kroku.R
# -----------------------------------------------------------------------------

deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi

rankine_krok_po_kroku <- function(
  q = 10.0,          # [kN/m2]
  gamma = 18.0,      # [kN/m3]
  phi_st = 30.0,     # [deg]
  epsilon_st = 0.0,  # [deg]
  gamma_Q = 1.0,     # [-]
  gamma_G = 1.0,     # [-]
  z_max = 5.0        # [m]
) {
  cat("=== Rankine: obliczenia krok po kroku ===\n\n")

  # 1) Dane i konwersja katow
  phi <- deg2rad(phi_st)
  epsilon <- deg2rad(epsilon_st)
  cat("1) Dane wejsciowe:\n")
  cat(sprintf("   q = %.3f kN/m2\n", q))
  cat(sprintf("   gamma = %.3f kN/m3\n", gamma))
  cat(sprintf("   phi = %.3f deg\n", phi_st))
  cat(sprintf("   epsilon = %.3f deg\n", epsilon_st))
  cat(sprintf("   gamma_Q = %.3f\n", gamma_Q))
  cat(sprintf("   gamma_G = %.3f\n", gamma_G))
  cat(sprintf("   z_max = %.3f m\n\n", z_max))

  # 2) Wartosci obliczeniowe
  q_d <- q * gamma_Q
  gamma_d <- gamma * gamma_G
  cat("2) Wartosci obliczeniowe:\n")
  cat(sprintf("   q_d = q * gamma_Q = %.6f kN/m2\n", q_d))
  cat(sprintf("   gamma_d = gamma * gamma_G = %.6f kN/m3\n\n", gamma_d))

  # 3) Kat omega_epsilon
  sin_ratio <- sin(epsilon) / sin(phi)
  if (abs(sin_ratio) > 1) {
    stop("Brak rozwiazania: |sin(epsilon)/sin(phi)| > 1. Sprawdz katy.")
  }
  omega_epsilon <- asin(sin_ratio)
  cat("3) Kat pomocniczy:\n")
  cat(sprintf("   omega_epsilon = asin(sin(epsilon)/sin(phi)) = %.6f rad = %.3f deg\n\n",
              omega_epsilon, rad2deg(omega_epsilon)))

  # 4) K_a oraz h_z_d (jak w kodzie Python)
  if (abs(epsilon) < 1e-6) {
    K_a <- tan(pi / 4 - phi / 2)^2
    h_z_d <- q_d / gamma_d
    cat("4) Dla epsilon ~ 0:\n")
    cat("   K_a = tan(pi/4 - phi/2)^2\n")
    cat("   h_z_d = q_d / gamma_d\n")
  } else {
    K_a <- sin(omega_epsilon - epsilon) * cos(epsilon) / sin(omega_epsilon + epsilon)
    h_z_d <- q_d / (gamma_d * cos(epsilon))
    cat("4) Dla epsilon != 0:\n")
    cat("   K_a = sin(omega_epsilon - epsilon) * cos(epsilon) / sin(omega_epsilon + epsilon)\n")
    cat("   h_z_d = q_d / (gamma_d * cos(epsilon))\n")
  }
  cat(sprintf("   K_a = %.6f\n", K_a))
  cat(sprintf("   h_z_d = %.6f m\n\n", h_z_d))

  # 5) Kata theta i 90 - theta
  theta <- pi / 4 + phi / 2 + (omega_epsilon - epsilon) / 2
  cat("5) Geometria klina:\n")
  cat(sprintf("   theta = %.6f rad = %.3f deg\n", theta, rad2deg(theta)))
  cat(sprintf("   90 - theta = %.3f deg\n\n", 90 - rad2deg(theta)))

  # 6) Rozklad parcia Rankine'a
  E_a <- K_a * gamma_d * (z_max^2 / 2 + h_z_d * z_max)
  E_ah <- E_a * cos(epsilon)
  E_av <- E_a * sin(epsilon)
  moment_1 <- K_a * gamma_d * (z_max^3 / 3 + h_z_d * z_max^2 / 2)
  z_c <- moment_1 / E_a

  cat("6) Sily i punkt przylozenia:\n")
  cat(sprintf("   E_a  = %.6f kN/m\n", E_a))
  cat(sprintf("   E_ah = %.6f kN/m\n", E_ah))
  cat(sprintf("   E_av = %.6f kN/m\n", E_av))
  cat(sprintf("   z_c  = %.6f m (od gory)\n\n", z_c))

  # 7) e_a(z) dla kilku glebokosci
  e_a <- function(z) K_a * (z + h_z_d) * gamma_d
  z_probe <- c(0, z_max / 2, z_max)
  cat("7) Przykladowe naprezenia jednostkowe e_a(z):\n")
  for (z in z_probe) {
    cat(sprintf("   e_a(%.3f m) = %.6f kN/m2\n", z, e_a(z)))
  }
  cat("\n=== Koniec obliczen ===\n")

  invisible(list(
    K_a = K_a,
    h_z_d = h_z_d,
    omega_epsilon_deg = rad2deg(omega_epsilon),
    theta_deg = rad2deg(theta),
    E_a = E_a,
    E_ah = E_ah,
    E_av = E_av,
    z_c = z_c
  ))
}

# Domyslne uruchomienie skryptu
wyniki <- rankine_krok_po_kroku()
