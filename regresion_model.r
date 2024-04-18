datos <- read.table("micro_milagro.csv", sep = ",", header = TRUE)

modelo <- lm(pes ~ inev, data=datos)

coefficients(modelo)

summary(modelo)

diagnorm(modelo)

