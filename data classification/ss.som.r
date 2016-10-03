## Self-Organizing Map (SOM) ##
install.packages("som")
library(som)
data.s = t(scale(t(data))) # row-wise scale : center to zero
data.som = som(data.s, xdim=4, ydim=4, topol="rect", neigh="gaussian")
	# Performs SOM clustering
	# topol = "hexa", "rect"

plot(data.som)
#data.som.group = data.som$visual$x*10+data.som$visual$y+1