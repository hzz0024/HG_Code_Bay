pop <- factor(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)), levels
              = c("Navarra", "Aragon", "Catalonia")) # Population
wing <- c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9) # Wingspan
body <- c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2) # Body length
sex <- factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F"))
mites <- c(0, 3, 2, 1, 0, 7, 0, 9, 6) # Number of ectoparasites
color <- c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57) # Color intensity
damage <- c(0,2,0,0,4,2,1,0,1) # Number of wings damaged
cbind(pop, sex, wing, body, mites, color, damage) # Print out data set

summary(fm1 <- lm(wing ~ pop + body))

summary(fm2 <- lm(wing ~ pop-1 + body))
fm2

coef(fm2)[4]

par(mfrow = c(1, 3), mar = c(5,4,2,2), cex = 1.2, cex.main = 1) # Figure 3.3 (a)
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14),
     lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Wingspan")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm2)[1], coef(fm2)[4], col = "red", lwd = 2)
abline(coef(fm2)[2], coef(fm2)[4], col = "blue", lwd = 2)
abline(coef(fm2)[3], coef(fm2)[4], col = "green", lwd = 2)
text(6.8, 14, "A", cex = 1.5)

model.matrix(~ pop*body)
model.matrix(~ pop*body-1)
model.matrix(~ pop*body-1-body)

summary(fm3 <- lm(wing ~ pop*body-1-body))
lm(formula = wing ~ pop * body - 1 - body)
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14),
     lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm3)[1], coef(fm3)[4], col = "red", lwd = 2)
abline(coef(fm3)[2], coef(fm3)[5], col = "blue", lwd = 2)
abline(coef(fm3)[3], coef(fm3)[6], col = "green", lwd = 2)
text(6.8, 14, "B", cex = 1.5)


model.matrix(~ pop+sex)
model.matrix(~ pop+sex-1)
summary(fm6 <- lm(wing ~ pop + sex-1))

model.matrix(~ pop*sex)
model.matrix(~ pop + sex + pop:sex)

model.matrix(~ body + color)
summary(fm7 <- lm(wing ~ body + color))


model.matrix(~ body*color)

summary(fm8 <- lm(wing ~ body + color + body*color))
summary(fm8 <- lm(wing ~ body*color))

body2 <- body^2
body3 <- body^3
model.matrix(~ body + body2 + body3)
summary(fm9 <- lm(wing ~ body + body2 + body3))
summary(fm9 <- lm(wing ~ body + I(body^2) + I(body^3))) # same

summary(fm10 <- glm(mites ~ pop-1 + body, family = poisson))

presence <- ifelse(mites > 0, 1, 0)
summary(fm11 <- glm(presence ~ pop-1 + body, family = binomial))


summary(fm12 <- glm(cbind(damage, 4-damage) ~ pop + body -1, family = binomial))

