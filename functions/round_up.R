##############
# Tanya Major
# April 2017
# round_up.R
#

round.up <- function(x, decimals = 1){
  number <- x + (5 * 10 ^ (-decimals - 1))
  round(number, digits = decimals)
}

round.down <- function(x, decimals = 1){
  number <- x - (5 * 10 ^ (-decimals - 1))
  round(number, digits = decimals)
}


