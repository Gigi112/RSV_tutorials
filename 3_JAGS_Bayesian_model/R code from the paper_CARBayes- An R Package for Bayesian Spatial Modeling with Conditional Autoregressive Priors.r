################################################################
#### Load the packages required to do an extended analysis
################################################################
library("foreign")
library("shapefiles")
library("sp")
library("boot")
library("Matrix")
library("nlme")
library("maptools")
library("grid")
library("deldir")
library("splines")
library("spdep")
library("CARBayes")
set.seed(1)


#############################
#### Example 1 - house prices
#############################
#### Read in the data
data("housedata", package = "CARBayes")
data("shp", package = "CARBayes")
data("dbf", package = "CARBayes")


#### Remove the outlying observation
housedata <- housedata[!rownames(housedata) == "S02000655", ]

#### Combine the data and shapefile
data.combined <- combine.data.shapefile(housedata, shp, dbf)

#### Plot a map of median property price
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000, 
  647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000, 
  647000), scale = 10000, fill = c("transparent", "black"))
text1 <- list("sp.text", c(225000, 649000), "0")
text2 <- list("sp.text", c(230000, 649000), "5000 m")
spplot(data.combined, "price", sp.layout = list(northarrow, scalebar, text1, text2), 
  scales = list(draw = TRUE), at = seq(min(housedata$price) - 1, max(housedata$price) + 
    1, length.out = 8), col.regions = c("#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", 
    "#EF3B2C", "#CB181D", "#99000D"))



#### Transform the price and crime variables
housedata$logprice <- log(housedata$price)
housedata$logcrime <- log(housedata$crime)
housedata$logdriveshop <- log(housedata$driveshop)

#### Fit a covariate only model
form <- "logprice ~ logcrime + rooms + sales + factor(type) + logdriveshop"
model <- lm(formula = form, data = housedata)
summary(model)

#### Compute a Moran's I test on the residuals
W.nb <- poly2nb(data.combined, row.names = rownames(housedata))
W.list <- nb2listw(W.nb, style = "B")
resid.model <- residuals(model)
moran.mc(x = resid.model, listw = W.list, nsim = 10000)


#### Run a regression model with spatially correlated random effects
W.mat <- nb2mat(W.nb, style = "B")
model.spatial <- gaussian.properCAR(as.formula(form), data = housedata,
  W = W.mat, burnin = 20000, n.sample = 1e+05, thin = 10)
model.spatial
summary(model.spatial)

#### Plot the spatial autocorrelation parameter
plot(model.spatial$samples$rho)




################################
#### Example 2 - disease mapping
################################
#### Read in the data
data("respdata", package = "CARBayes")
data("shp", package = "CARBayes")
data("dbf", package = "CARBayes")


#### Combine the data and shapefile and create the neighborhood matrix W
respdata$SIR2010 <- respdata$observed2010/respdata$expected2010
data.combined <- combine.data.shapefile(respdata, shp, dbf)
W.nb <- poly2nb(data.combined, row.names = rownames(respdata))
W.mat <- nb2mat(W.nb, style = "B")


#### Plot a map of the SIR
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000, 
  647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000, 
  647000), scale = 10000, fill = c("transparent", "black"))
text1 <- list("sp.text", c(225000, 649000), "0")
text2 <- list("sp.text", c(230000, 649000), "5000 m")
spplot(data.combined, "SIR2010", sp.layout = list(northarrow, scalebar, text1, 
  text2), scales = list(draw = TRUE), at = seq(min(respdata$SIR2010) - 0.05, max(respdata$SIR2010) + 
  0.05, length.out = 8), col.regions = c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", 
  "#FC4E2A", "#E31A1C", "#B10026"))



#### Create the dissimilarity metric
Z.incomedep <- as.matrix(dist(cbind(respdata$incomedep2010, respdata$incomedep2010), 
  method = "manhattan", diag = TRUE, upper = TRUE)) * W.mat/2

#### Run the local CAR model
formula <- "observed2010 ~ offset(log(expected2010))"
model.dissimilarity <- poisson.dissimilarityCAR(as.formula(formula), data = respdata,
  W = W.mat, Z = list(Z.incomedep = Z.incomedep), 
  rho = 0.99, fix.rho = TRUE, burnin = 20000, n.sample = 1e+05, thin = 10)
model.dissimilarity

#### Plot a map with the boundaries overlaid
border.locations <- model.dissimilarity$W.summary$W.posterior
risk.estimates <- model.dissimilarity$fitted.values[, 3]/respdata$expected2010
data.combined@data <- data.frame(data.combined@data, risk.estimates)
boundary.final <- highlight.borders(border.locations = border.locations, ID = rownames(respdata), 
  shp = shp, dbf = dbf)
boundaries <- list("sp.points", boundary.final, col = "white", pch = 19, cex = 0.2)
northarrow <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000, 
  647000), scale = 4000)
scalebar <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000, 
  647000), scale = 10000, fill = c("transparent", "black"))
text1 <- list("sp.text", c(225000, 649000), "0")
text2 <- list("sp.text", c(230000, 649000), "5000 m")
spplot(data.combined, "risk.estimates", sp.layout = list(northarrow, scalebar, 
  text1, text2, boundaries), scales = list(draw = TRUE), at = seq(min(risk.estimates) - 
  0.1, max(risk.estimates) + 0.1, length.out = 8), col.regions = c("#FFFFB2", "#FED976", 
  "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026")) 
