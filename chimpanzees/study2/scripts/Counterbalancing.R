library("readxl")

left.location = c("high", "low", "none")
middle.location = c("high", "low", "none")
right.location = c("high", "low", "none")



# create balanced predictors
counterbalancing <- data.frame(expand.grid(
  left.location = left.location,
  middle.location = middle.location, 
  right.location = right.location,
  available = c("left - middle", "middle - right", "left - right")
))


write.table(counterbalancing, file = "./counterbalancing.txt", sep="\t", row.names = TRUE)
