library(sp)
library(sf)
datos <- ngs.aea
datos <- spTransform(ngs.aea,CRS('+proj=laea +lat_0=90 +lon_0=-150 +x_0=0 +
                                 y_0=0 +datum=WGS84 +units=m +no_defs'))

alaska<- subset(datos, CATEGORY == 'SW-ALASKA')
alaska<- subset(alaska, ZN_ICP40 > 0)
str(alaska)

zero <- zerodist(alaska)
length(unique(zero[,1]))
alaska <- remove.duplicates(alaska) #se eliminan duplicados

# define groups for mapping
quantile(alaska$ZN_ICP40, prob=c(0,0.25,0.5,0.75,1))
cuts <-c(0,64,78,95,418,450)
# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'darkorange', 'blue', 'dark blue'))
spplot(alaska, 'ZN_ICP40', cuts=cuts, col.regions=blues(5),  pch=20, cex=2,
       main="Suroeste de Alaska: Concentración del zinc por quintiles")

#se ve la distribución de la variable
##Variograma
library(gstat)
zn <- gstat(formula=ZN_ICP40~1, locations=alaska)
n <- variogram(zn)
summary(n)
plot(n, main="Semivariograma de la concetración de zinc",
     xlab="Distancia",ylab="Semivarianza",
     pch=20,
     cex=1.5)
plot(variogram(ZN_ICP40~1, alaska, alpha = c(0, 45, 90, 135)),
     pch=20,
     cex=1.5,
     xlab="Distancia", ylab="Semivarianza",
     main="Semivariograma de la concetración de zinc en dirección a los puntos cardinales")
#hacer el variograma en varias direcciones
###
par(mfrow=c(2,2))
fS <- fit.variogram(n, vgm("Sph"), fit.kappa = TRUE)
plot(variogramLine(fS, 300000), type='l',ylim=c(0,1200) ,col='blue', lwd=2, main = "Esférico",xlab="Distancia", ylab="Semivarianza")
points(n[,2:3], pch=20, col='red')

fE <- fit.variogram(n, vgm("Exp"), fit.kappa = TRUE)
plot(variogramLine(fE, 300000), type='l',ylim=c(0,1200) ,col='blue',lwd=2, main = "Exponencial",xlab="Distancia", ylab="Semivarianza")
points(n[,2:3], pch=20, col='red')

fG <- fit.variogram(n, vgm("Gau"), fit.kappa = TRUE)
plot(variogramLine(fG, 300000), type='l',ylim=c(0,1200) ,col='blue', lwd=2, main = "Gaussiano",xlab="Distancia", ylab="Semivarianza")
points(n[,2:3], pch=20, col='red')

fM <- fit.variogram(n, vgm("Mat"), fit.kappa = TRUE)
plot(variogramLine(fM, 300000), type='l',ylim=c(0,1200) ,col='blue', lwd=2, main = "Matern",xlab="Distancia", ylab="Semivarianza")
points(n[,2:3], pch=20, col='red')
par(mfrow=c(1,1))
#Encontrar el modelo que mejor ajusta
fit.variogram(n, vgm(c("Exp", "Mat", "Sph","Gau")), fit.kappa = TRUE)
# es el modelo Matern
library(raster)
#Kriging con gstat
r <- raster(alaska, nrows=300, ncols=200)
res(r)
g <- as(r, 'SpatialGrid')

ke_E <- krige(formula = ZN_ICP40 ~ 1 ,locations = alaska, newdata=g, model=fE)
ke_S <- krige(formula = ZN_ICP40 ~ 1 ,locations = alaska, newdata=g, model=fS)
ke_G <- krige(formula = ZN_ICP40 ~ 1 ,locations = alaska, newdata=g, model=fG)
ke_M <- krige(formula = ZN_ICP40 ~ 1 ,locations = alaska, newdata=g, model=fM)

plot(ke_E, main= "Exponencial")
plot(ke_S, main = "Esférico")
plot(ke_G, main = "Gaussiano")
plot(ke_M, main = "Matern")

