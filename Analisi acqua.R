# Analisi acqua

# UPLOAD
library(geoR)
temp_acqua<-read.table("temp_acqua.dat")
is.geodata(temp_acqua)
geotemp <- as.geodata(temp_acqua,coords=4:5,data.col=3) # perché le coordinate sono in 4 e 5 colonna


# ANALISI ESPLORATIVA SPAZIALE, con commento dei risultati
head(temp_acqua)
summary(temp_acqua)
names(temp_acqua)

# eseguiamo il grafico a cerchi
points(geotemp, cex.min = 0.2, cex.max = 1, pt.divide = 4)
points(geotemp, cex.min = 0.2, cex.max = 1, pt.divide = 4, trend = "1st")
points(geotemp, cex.min = 0.2, cex.max = 1, pt.divide = 4, trend = "2nd")

# c'è una dipendenza dallo spazio, poiché i cerchi piu grandi sono collocati piu a sud

plot(geotemp, lowess = TRUE)
# sono presenti due possibili trend per le due coordinate, inoltre, rispetto alla curva di lisciamento,
# per la coordinata X rispetto ai dati c'è un punto che si discosta molto dalla massa,
# tuttavia appunto per entrambe le coordinate sembrano esserci delle componenti di trend
# il grafico della distribuzione empirica dei dati risulta asimmetrico, tuttavia un'eventuale eliminazione delle osservazioni
# nella coda di destra potrebbe risultare in una distribuzione più campanulare, per tanto l'applicazione futura
# di un modello gaussiano non dovrebbe essere scartata a priori


# VARIOGRAMMA CAMPIONARIO
v0 = variog(geotemp)
plot(v0, type = "b", lwd = 2) # lanciamo il variogramma

# il testo ora ci dice di considerare una distanza massima pari a 50 con un passo massimo di 5 km, e noi lo facciamo:
v00 = variog(geotemp, uvec = seq(0,50, by=5))
plot(v00, type = "b", lwd = 2, main = "Variogramma campionario")
# 50 risulta una distanza sensata, perché per quel valore il variogramma raggiunge il sill (noto che per valori
# maggiori, ad esempio 80, vi è una repentina decrescita e ri-crescita della correlazione spaziale proprio dopo quel valore, quindi
# è lecito presumere che il sill venga raggiunto proprio in quel punto)
# una restrizione dell'intervallo rende più difficile la lettura, mettendo 7 invece
# acquisiamo più informazione dalla sua lettura
v00 = variog(geotemp, uvec = seq(0,50,by=7))
plot(v00, type = "b", lwd = 2, main = "Variogramma campionario")


# TREND PER IL VARIOGRAMMA CAMPIONARIO
# Costruiamo ora un variogramma sui RESIDUI di un modello per TREND LINEARE della media
# perché lineare? perché c'era un'evidenza di linearità nello scatterplot tra coordinata Y e dati e X e dati
# SE RISPETTO ALLE COORDINATE LA NUVOLA DI PUNTI FOSSE DISPERSA OMOGENEAMENTE, AVREMMO ASSUNTO TREND COSTANTE, ma
# non è questa la situazione
v_lin_res = variog(geotemp, uvec = seq(0, 50, by = 5), trend = "1st")
plot(v_lin_res, type = "b", lwd = 2)

# anziché plottarlo e basta, ovviamente, essendo che sono definiti sulla stessa scala,
# lo sovrapponiamo al grafico precedente e vediamo quale è migliore in termini di livello di 
# correlazione spaziale
lines(v_lin_res, type = "b", lwd = 2, col = "red")
# il variogramma campionario basato sui residui di un trend lineare per la media, è migliore
# in quanto la correlazione spaziale risulta diminuita, MA QUI IL VARIOGRAMMA SI STABILIZZA A UN RANGE
# PARI A 12 CIRCA! Tra 0 e 50 avevamo notato altresì che il modello non risulta
# stazionario, questo è importante perché la stazionarietà è collegata alla bontà delle stime e dell'adattamento
# questo ci dice, secondo il modello per il trend nella media, che se noi all'interno della regione Veneto
# ci spostiamo per distanze fino a un massimo di 12 km, la temperatura risulta fortemente correlata

v_quad_res = variog(geotemp, uvec = seq(0, 50, by = 7), trend = "2nd")
lines(v_quad_res, type = "b", lwd = 2, col = "blue")
# il modello per trend quadratico viene escluso in quanto, seppur non di molto, aumenta la correlazione spaziale
# per distanze pari a 12km, lasciando quindi intuire che il variogramma basato sui residui di un trend lineare per
# la media sia il migliore per tentare di descrivere il fenomeno


# VARIOGRAMMA PARAMETRICO
# avendo visto che il modello per trend lineare è il migliore, studiamo quale famiglia parametrica
# si adatta meglio a QUESTO specifico modello
variofit_mat = variofit(v_lin_res, ini = c(5, 25), cov.model = "matern", fix.nugget = TRUE)
variofit_gau = variofit(v_lin_res, ini = c(5, 25), cov.model = "gaussian", fix.nugget = TRUE)
variofit_sph = variofit(v_lin_res, ini = c(5, 25), cov.model = "spheric", fix.nugget = TRUE)
# non ho specificato i valori iniziali! Il software ha provato quindi a farlo da solo (il famoso 'ini'...), 
# questi dipendono dal variogramma campionario ovviamente

# plottiamo ora il variogramma campionario che avevamo scelto, ovvero quello a trend
# lineare, e gli sovrapponiamo il modello con famiglia Matern, in modo da studiarne l'adattamento:
plot(v_lin_res, lwd = 2, type = "b")
lines.variomodel(variofit_mat, col = "red", lwd = 2)

# aggiungiamo quello gaussiano:
lines.variomodel(variofit_gau, col = "blue", lwd = 2)

# e infine quello sferico:
lines.variomodel(variofit_sph, col = "dark green", lwd = 2)
# visivamente, l'adattamento migliore � quello del modello di Matern, in quanto l'interpolazione
# al variogramma campionario risulta migliore degli altri modelli
# solo il modello matern tra l'altro ha effetto nugget nullo, nel modello gaussiano e 
# sferico c'è effetto nugget

variofit_mat05 = variofit(v_lin_res, ini = c(5, 25), cov.model = "matern", fix.nugget = TRUE, kappa = 0.5)
variofit_mat15 = variofit(v_lin_res, ini = c(5, 25), cov.model = "matern", fix.nugget = TRUE, kappa = 1.5)
variofit_mat25 = variofit(v_lin_res, ini = c(5, 25), cov.model = "matern", fix.nugget = TRUE, kappa = 2.5)
plot(v_lin_res, lwd = 2, type = "b")
lines.variomodel(variofit_mat05, col = "red", lwd = 2)
lines.variomodel(variofit_mat15, col = "blue", lwd = 2)
lines.variomodel(variofit_mat25, col = "orange", lwd = 2)
variofit_mat05
variofit_mat15
# il modello 1.5 ha minor varianza e minor range, lo preferiamo

# MODELLO GEOSTATISTICO GAUSSIANO
hist(temp_acqua$Temp)
# 0.5,  1,  1.5  sulla base del bic
mlgau05 = likfit(geotemp, fix.nugget = TRUE, ini = c(5, 25), cov.model = "matern", kappa = 0.5, trend = "1st")
mlgau05
mlgau1 = likfit(geotemp, fix.nugget = TRUE, ini = c(5, 25), cov.model = "matern", kappa = 1, trend = "1st")
mlgau1
mlgau15 = likfit(geotemp, fix.nugget = TRUE, ini = c(5, 25), cov.model = "matern", kappa = 1.5, trend = "1st")
mlgau15
plot(v_lin_res, lwd = 2, type = "b")
lines.variomodel(seq(0, 50, l = 1000), cov.pars = mlgau05$cov.pars, cov.model = "matern", kappa = 0.5, nugget = mlgau05$nugget, col = "red", lwd = 2)
lines.variomodel(seq(0, 50, l = 1000), cov.pars = mlgau1$cov.pars, cov.model = "matern", kappa = 1, nugget = mlgau1$nugget, col = "blue", lwd = 2)
lines.variomodel(seq(0, 50, l = 1000), cov.pars = mlgau15$cov.pars, cov.model = "matern", kappa = 1.5, nugget = mlgau15$nugget, col = "orange", lwd = 2)
legend("bottomright", c("0.5","1","1.5"), col = c("red", "blue", "orange"))

# kriging stima sul punto Padova x = 725, y = 5032
plot(geotemp)
locs = pred_grid(c(600, 800), c(5000, 5200), by = 2)
KC = krige.control(type = "SK", obj.model = mlgau15)
sk = krige.conv(geotemp, krige = KC, locations = locs)
sk$predict
sk$krige.var
pred.lim = range(sk$predict)
image(sk, col = gray(seq(1, 0, l = 50)), zlim = pred.lim, main = "kriging estimates")
contour(sk, add = TRUE, nlev = 6)
points(geotemp, add = TRUE, cex.max = 2, col = "red")