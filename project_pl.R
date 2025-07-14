#===== 1. Wczytanie danych =====
library(readxl)
library(ggplot2)
library(numDeriv)

Dane <- read_excel("Dane 8.1.xlsx")


#===== Do generowania danych =====
# n <- 200  # liczba na grupę
# 
# alpha0 <- 1.5
# lambda0 <- 8
# 
# alpha1 <- 2
# lambda1 <- 6
# 
# med0 <- lambda0 * (log(2))^(1/alpha0)
# med1 <- lambda1 * (log(2))^(1/alpha1)
# cat("Teoretyczna mediana G=0:", round(med0, 2), "\n")
# cat("Teoretyczna mediana G=1:", round(med1, 2), "\n")
# 
# T0 <- rweibull(n, shape = alpha0, scale = lambda0)
# T1 <- rweibull(n, shape = alpha1, scale = lambda1)
# 
# # Cenzorowanie: np. 20% przypadków
# D0 <- rbinom(n, 1, 0.8)
# D1 <- rbinom(n, 1, 0.8)
# 
# # Składamy dane
# T_all <- c(T0, T1)
# D_all <- c(D0, D1)
# G_all <- c(rep(0, n), rep(1, n))
# 
# Dane <- data.frame(T = T_all, D = D_all, G = G_all)



#===== 2. Przygotowanie danych =====
T <- Dane$T
sum(T==0) # 2
T[T == 0] <- 0.0025  # było kilka zer
sum(T==0)

D <- Dane$D
G <- Dane$G



#===== 3. Funkcja do wyznaczenia parametrów Weibulla =====
lik_dane <- function(x) {
  alpha1 <- x[1]  # shape dla G = 1
  lambda1 <- x[2] # scale dla G = 1
  alpha0 <- x[3]  # shape dla G = 0
  lambda0 <- x[4] # scale dla G = 0
  
  grupa1 <- (G == 1) * (D * (log(alpha1 / lambda1) + (alpha1 - 1) * log(T / lambda1)) - (T / lambda1)^alpha1)
  grupa0 <- (G == 0) * (D * (log(alpha0 / lambda0) + (alpha0 - 1) * log(T / lambda0)) - (T / lambda0)^alpha0)
  
  return(-sum(grupa1 + grupa0))  
}



#===== 4. Estymacja parametrów czyli MLE =====
mle_dane <- optim(par = c(1, 1, 1, 1), fn = lik_dane, method = "L-BFGS-B", lower = rep(1e-6, 4))

# method = "L-BFGS-B" i lower = rep(1e-6, 4) 
# aby zapobiec ujemnym wartościom parametrów



#===== 5. Wyodrębnienie oszacowanych parametrów =====
alpha1 <- mle_dane$par[1]  # ≈ 0.8309
lambda1 <- mle_dane$par[2] # ≈ 4.6659
alpha0 <- mle_dane$par[3]  # ≈ 1.8774
lambda0 <- mle_dane$par[4] # ≈ 3.8744



#===== 6. Do wykresu =====
x <- seq(0, 30, length.out = 300)

df_plot <- data.frame(
  t = rep(x, 2),
  S = c(pweibull(x, shape = alpha1, scale = lambda1, lower.tail = FALSE),
        pweibull(x, shape = alpha0, scale = lambda0, lower.tail = FALSE)),
  grupa = factor(rep(c("G = 1", "G = 0"), each = length(x)))
)



#===== 7. Wykres z ggplot2 =====
ggplot(df_plot, aes(x = t, y = S, color = grupa)) +
  geom_line(size = 1.2) +
  labs(title = "Krzywe przeżycia Weibulla",
       x = "Czas",
       y = "Prawdopodobieństwo przeżycia",
       color = "Grupa") +
  theme_minimal() +
  theme(text = element_text(size = 13))


#===== 8. Obliczanie median przeżycia dla obu grup =====
mediana_G1 <- lambda1 * (log(2))^(1 / alpha1)
mediana_G0 <- lambda0 * (log(2))^(1 / alpha0)

cat("Mediana G = 1:", round(mediana_G1, 3),
    "\nMediana G = 0:", round(mediana_G0, 3))



#===== 9. Macierz wariancji (asymptotyczna) =====
# Funkcja lik() jest ujemną log-wiarygodnością
# Więc macierz Hessego = - macierz informacji Fishera (obserwowaneJ)
# library(numDeriv)
hessian_macierz <- hessian(func = lik_dane, x = mle_dane$par)
vcov_matrix <- solve(hessian_macierz)
# Otrzymujemy asymptotyczną macierz wariancji-kowariancji estymatorów 



#===== 10. Metoda delta – PU dla G = 1 =====
# Interesujący nas fragment macierzy kowariancji
V1 <- vcov_matrix[1:2, 1:2]

# Pochodne do gradientu
# mediana = g(alpha, lambda) = lambda * (log(2))^(1/alpha)
d_alpha1 <- -lambda1 * (log(2))^(1 / alpha1) * log(log(2)) / (alpha1^2)
d_lambda1 <- (log(2))^(1 / alpha1)

grad1 <- c(d_alpha1, d_lambda1)

# Wariancja = gradient^T * V * gradient
var_mediana1 <- t(grad1) %*% V1 %*% grad1
se_mediana1 <- sqrt(var_mediana1)

# Przedziały ufności
PU1 <- c(mediana_G1 - 1.96 * se_mediana1, mediana_G1 + 1.96 * se_mediana1)



#===== 11. Metoda delta – PU dla G = 0 =====
# To samo co przedtem ale dla drugich parametrów
V0 <- vcov_matrix[3:4, 3:4]
d_alpha0 <- -lambda0 * (log(2))^(1 / alpha0) * log(log(2)) / (alpha0^2)
d_lambda0 <- (log(2))^(1 / alpha0)
grad0 <- c(d_alpha0, d_lambda0)
var_mediana0 <- t(grad0) %*% V0 %*% grad0
se_mediana0 <- sqrt(var_mediana0)
PU0 <- c(mediana_G0 - 1.96 * se_mediana0, mediana_G0 + 1.96 * se_mediana0)



#===== 12. Otrzymane przedziały ufności 95% =====
cat("Mediana G = 1:", round(mediana_G1, 3),
    "\n95% PU dla mediany G = 1:", round(PU1[1], 3), "-", round(PU1[2], 3), "\n",
    "\nMediana G = 0:", round(mediana_G0, 3),
    "\n95% PU dla mediany G = 0:", round(PU0[1], 3), "-", round(PU0[2], 3))



#===== 13. Wykres z zaznaczonymi medianami i przedziałami ufności=====
x <- seq(2, 4, length.out = 300) # gdy nne dane - zmienić końce przedzuału !!!

df_plot <- data.frame(
  t = rep(x, 2),
  S = c(pweibull(x, shape = alpha1, scale = lambda1, lower.tail = FALSE),
        pweibull(x, shape = alpha0, scale = lambda0, lower.tail = FALSE)),
  grupa = factor(rep(c("G = 1", "G = 0"), each = length(x))))

ggplot(df_plot, aes(x = t, y = S, color = grupa)) +
  geom_line(size = 1.2) +
  # Pozioma linia S = 0.5
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  # Linie pionowe – mediany
  geom_vline(xintercept = mediana_G1, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = mediana_G0, linetype = "dashed", color = "red") +
  # Linie pionowe – PU dla mediany (95%)
  geom_vline(xintercept = PU1[1], linetype = "dotted", color = "blue") +
  geom_vline(xintercept = PU1[2], linetype = "dotted", color = "blue") +
  geom_vline(xintercept = PU0[1], linetype = "dotted", color = "red") +
  geom_vline(xintercept = PU0[2], linetype = "dotted", color = "red") +
  labs(title = "Krzywe przeżycia Weibulla z medianami i 95% PU",
       x = "Czas",
       y = "Prawdopodobieństwo przeżycia",
       color = "Grupa") +
  theme_minimal() +
  theme(text = element_text(size = 13))
  


#===== 14. Bootstrap dla mediany w obu grupach =====
# Inny sposób obliczenia przedziałów ufności

B <- 1000  # liczba powtórzeń

# liczba generowanych próbek dla obu grup
n1 <- sum(G == 1)
n0 <- sum(G == 0)

# tablice do zapisania wyników
boot_median_G1 <- numeric(B)
boot_median_G0 <- numeric(B)

# Bootstrap G = 1
for (b in 1:B) {
  T_sim <- rweibull(n1, shape = alpha1, scale = lambda1)
  # MLE dla danych symulowanych
  mle <- suppressWarnings(
    optim(par = c(1, 1), 
    fn = function(x) {
      a <- x[1] 
      l <- x[2]
      -sum(log(a / l) + (a - 1) * log(T_sim / l) - (T_sim / l)^a)
    }, 
    method = "L-BFGS-B", lower = rep(1e-6, 2))
  )
  a_hat <- mle$par[1]
  l_hat <- mle$par[2]
  boot_median_G1[b] <- l_hat * (log(2))^(1 / a_hat)
}

# Bootstrap G = 0
for (b in 1:B) {
  T_sim <- rweibull(n0, shape = alpha0, scale = lambda0)
  mle <- suppressWarnings(
    optim(par = c(1, 1), 
    fn = function(x) {
      a <- x[1]
      l <- x[2]
      -sum(log(a / l) + (a - 1) * log(T_sim / l) - (T_sim / l)^a)
    }, 
    method = "L-BFGS-B", lower = rep(1e-6, 2))
  )
  a_hat <- mle$par[1]
  l_hat <- mle$par[2]
  boot_median_G0[b] <- l_hat * (log(2))^(1 / a_hat)
}

# Przedziały ufności (percentylowe)
PU_boot_G1 <- quantile(boot_median_G1, probs = c(0.025, 0.975))
PU_boot_G0 <- quantile(boot_median_G0, probs = c(0.025, 0.975))

# Wyniki
cat("Mediana G = 1:", round(mediana_G1, 3),
    "\n95% PU dla mediany G = 1:", round(PU1[1], 3), "-", round(PU1[2], 3),
    "\nBootstrap PU dla mediany G = 1:", round(PU_boot_G1, 3), "\n",
    "\nMediana G = 0:", round(mediana_G0, 3),
    "\n95% PU dla mediany G = 0:", round(PU0[1], 3), "-", round(PU0[2], 3),
    "\nBootstrap PU dla mediany G = 0:", round(PU_boot_G0, 3), "\n")



#===== 15. Test ilorazu wiarygodności =====

# (H1) – bez ograniczenia mediany - było już

mle_H1 <- optim(par = c(1, 1, 1, 1), fn = lik_dane, method = "L-BFGS-B", lower = rep(1e-6, 4))
mle_H1$par
lrtest_H1 <- mle_H1$value

# Dopasowanie modelu z narzuconą wspólną medianą (H0)
lik_H0 <- function(x) {
  alpha1 <- x[1]
  lambda1 <- x[2]
  alpha0 <- x[3]
  
  lambda0 <- lambda1 * (log(2))^(1/alpha1 - 1/alpha0)
  
  grupa1 <- (G == 1) * (D * (log(alpha1 / lambda1) + (alpha1 - 1) * log(T / lambda1)) - (T / lambda1)^alpha1)
  grupa0 <- (G == 0) * (D * (log(alpha0 / lambda0) + (alpha0 - 1) * log(T / lambda0)) - (T / lambda0)^alpha0)
  
  return(-sum(grupa1 + grupa0))
}

mle_H0 <- optim(par = c(1, 1, 1), fn = lik_H0, method = "L-BFGS-B", lower = rep(1e-6, 3))
mle_H0$par
lrtest_H0 <- mle_H0$value

# Statystyka LR
LR <- 2 * (lrtest_H0 - lrtest_H1)
pwartosc_lrtest <- 1 - pchisq(LR, df = 1) # 4-3 parametry

# 4. Wynik
cat("Statystyka LR:", round(LR, 3),
    "\np-wartość testu:", round(pwartosc_lrtest, 4))



#===== 16. Bootstrapowy test różnicy median pod H0: wspólna mediana =====
B <- 1000
n1 <- sum(G == 1)
n0 <- sum(G == 0)

# Dopasowanie Weibulla z narzuconą wspólną medianą
mle_boot_h0 <- optim(par = c(1, 1, 1), fn = function(x) {
  alpha1 <- x[1]  # shape dla G=1
  lambda1 <- x[2] # scale dla G=1
  alpha0 <- x[3]  # shape dla G=0
  
  # Wyliczamy scale dla G=0 tak, aby mediana była taka sama jak w G=1
  lambda0 <- lambda1 * (log(2))^(1 / alpha1 - 1 / alpha0)
  
  # Log-wiarygodność
  grupa1 <- (G == 1) * (D * (log(alpha1 / lambda1) + (alpha1 - 1) * log(T / lambda1)) - (T / lambda1)^alpha1)
  grupa0 <- (G == 0) * (D * (log(alpha0 / lambda0) + (alpha0 - 1) * log(T / lambda0)) - (T / lambda0)^alpha0)
  
  return(-sum(grupa1 + grupa0))
}, method = "L-BFGS-B", lower = rep(1e-6, 3))

alpha_boot_h0_1 <- mle_boot_h0$par[1]
lambda_boot_h0_1 <- mle_boot_h0$par[2]
alpha_boot_h0_0 <- mle_boot_h0$par[3]
lambda_boot_h0_0 <- lambda_boot_h0_1 * (log(2))^(1 / alpha_boot_h0_1 - 1 / alpha_boot_h0_0)

# Symulacja pod H0
boot_roznice_h0 <- numeric(B)

for (b in 1:B) {
  # Generujemy próbki z Weibulla z tą samą medianą, ale różnymi kształtami
  T_sim_G1 <- rweibull(n1, shape = alpha_boot_h0_1, scale = lambda_boot_h0_1)
  T_sim_G0 <- rweibull(n0, shape = alpha_boot_h0_0, scale = lambda_boot_h0_0)
  
  # Dopasowanie Weibulla do każdej symulowanej grupy
  mle1 <- optim(par = c(1, 1), fn = function(par) {
    a <- par[1]; l <- par[2]
    -sum(log(a / l) + (a - 1) * log(T_sim_G1 / l) - (T_sim_G1 / l)^a)
  }, method = "L-BFGS-B", lower = rep(1e-6, 2))
  
  mle0 <- optim(par = c(1, 1), fn = function(par) {
    a <- par[1]; l <- par[2]
    -sum(log(a / l) + (a - 1) * log(T_sim_G0 / l) - (T_sim_G0 / l)^a)
  }, method = "L-BFGS-B", lower = rep(1e-6, 2))
  
  m1 <- mle1$par[2] * (log(2))^(1 / mle1$par[1])
  m0 <- mle0$par[2] * (log(2))^(1 / mle0$par[1])
  
  boot_roznice_h0[b] <- m0 - m1
}

# Obserwowana różnica
obs_roznica <- mediana_G0 - mediana_G1

# p-wartość 
pwartosc_boot <- mean(abs(boot_roznice_h0) >= abs(obs_roznica))
sum(abs(boot_roznice_h0) >= abs(obs_roznica))

abs(obs_roznica)
head(abs(boot_roznice_h0))

# Wynik
cat("Obserwowana różnica median:", round(abs(obs_roznica), 3),
    "\np-wartość testu:", round(pwartosc_boot, 4), "\n")







#######################################################################
test_W = (abs(mediana_G0-mediana_G1))/(sqrt(var_mediana1+var_mediana0))
pwartosc_W = 2*(1 - pnorm(abs(test_W)))

cat("Obserwowana różnica median:", round(abs(obs_roznica), 3),
    "\np-wartość testu W:", round(pwartosc_W, 4), "\n")
