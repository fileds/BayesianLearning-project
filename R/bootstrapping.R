library(bnlearn)
rm(list = ls())

df <- gaussian.test[1:200,]

n_params <- ncol(df)
n_bootstraps <- 1000
pdags <- lapply(1:n_bootstraps, function(x) NULL)
ws <- lapply(1:n_bootstraps, function(x) NULL)
w_avg <- matrix(rep(0, n_params*n_params), n_params, n_params)
for (i in 1:n_bootstraps)
{
  idxs <- as.integer(sample(1:200, 200, replace = TRUE))
  df_i <- df[idxs, ]
  pdag <- pc.stable(df_i)
  pdags[[i]] <- pdag
  ws[[i]] <- amat(pdag)
  w_avg <- w_avg + amat(pdag)
}
w_avg <- w_avg / n_bootstraps
w_avg

uws <- unique(ws)
# Calculate the number of occurancies of each matrix
occurancies <- tabulate(match(ws, uws)) / n_bootstraps
occurancies
sort(occurancies)
most_likely <- uws[[4]]
second_most_likely <- uws[[2]]
most_likely - second_most_likely
# Doesnt map to uws
graphviz.plot(pdags[[4]])
graphviz.plot(pdags[[2]])
true <- model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
graphviz.plot((true))

all.equal(uws[[4]], uws[[4]])
which(all.equal(uws[[4]], ws))

i <- 1
while (TRUE)
{
  print(i)
  status <- all.equal(uws[[2]], ws[[i]])
  if (status == TRUE)
  {
    break
  }
  i <- i + 1
}
ws[[i]]
uws[[4]]
all.