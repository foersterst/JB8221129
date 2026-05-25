# R function to calculate taxonomic, phylogenetic or functional beta diversity, plus variation partitioning and significance tests
# Stenio I A Foerster
# https://foersterst.github.io/
# April 12, 2026

# comm: community matrix (sites x species). Abundances will be converted to presence-absence
# envs: environmental variables, can be a mix of continuous and discrete variables
# mems: distance-based Moran eigenvector maps
# tree: a phylogenetic (class phylo) tree or dendogram built from functional traits
# response: string; "total", "repl", "rich", "gain", or "loss"


fit_beta_div <- function(comm, envs, mems, tree, response) {
  # prepare input objects ----
  comm <- as.data.frame(comm)
  envs <- as.data.frame(envs)
  mems <<- as.data.frame(mems)

  # beta diversity object ----
  if (missing(tree)) {
    bt <- BAT::beta(comm = comm, func = "sorensen", abund = FALSE)
  } else {
    bt <- BAT::beta(comm = comm, tree = tree, func = "sorensen", abund = FALSE)
  }

  # summarize matrices ----
  if (response != "total") {
    ms <- "Summary available only for response = 'total'"
  } else {
    ms <- data.frame(
      "mean" = sapply(bt, mean, na.rm = TRUE),
      "min" = sapply(bt, min, na.rm = TRUE),
      "max" = sapply(bt, max, na.rm = TRUE)
    )
  }

  # response matrix ----
  if (response == "total") {
    y <- bt$Btotal
  } else if (response == "repl") {
    y <- bt$Brepl
  } else if (response == "rich") {
    y <- bt$Brich
  } else if (response == "gain") {
    y <- bt$Bgain
  } else if (response == "loss") {
    y <- bt$Bloss
  }

  # relative contributions ----
  x1 <- mean(bt$Btotal)
  x2 <- mean(bt$Brepl) / x1
  x3 <- mean(bt$Brich) / x1

  if (response == "total") {
    r_contrib <- c(x1, x2, x3)
  } else if (response == "repl") {
    r_contrib <- c(mean(bt$Brepl), NA, NA)
  } else if (response == "rich") {
    r_contrib <- c(mean(bt$Brich), NA, NA)
  } else if (response == "gain") {
    r_contrib <- c(mean(bt$Bgain), NA, NA)
  } else if (response == "loss") {
    r_contrib <- c(mean(bt$Bloss), NA, NA)
  }

  names(r_contrib) <- c("mean_matrix", "replacement", "rich_diff")

  # forward select MEMS ----
  gm <- dbrda(y ~ ., data = mems) # global model
  nm <- dbrda(y ~ 1, data = mems) # null model
  fw <- ordiR2step(object = nm, scope = formula(gm), permutations = 1000, trace = FALSE) # reduced model (selected MEMs)
  mems_string <- attr(fw$terms, "term.labels") # selected mems

  if (length(mems_string) > 0) {
    sel_mems <- mems[, mems_string]
  } else {
    mems_string <- c("No MEMs were selected. Proceeding with all MEMs.")
    sel_mems <- mems
  }

  # variation partitioning ----
  vp <- varpart(y, sel_mems, envs)

  # permutation tests ----
  xc <- data.frame(sel_mems, envs) # combine mems and envs
  perm_abc <- permutest(dbrda(y ~ ., data = xc), permutations = 1000) # fraction abc
  perm_a <- permutest(dbrda(y ~ ., data = sel_mems), permutations = 1000) # fraction a
  perm_b <- permutest(dbrda(y ~ ., data = envs), permutations = 1000) # fraction b
  # individual fractions
  a <- paste(colnames(sel_mems), collapse = " + ")
  b <- paste(colnames(envs), collapse = " + ")
  foA <- as.formula(paste("y ~", a, "+ Condition(", b, ")")) # formula for fraction A
  foB <- as.formula(paste("y ~", b, "+ Condition(", a, ")")) # formula for fraction B
  perm_pure_a <- permutest(dbrda(foA, data = xc), permutations = 1000) # pure a
  perm_pure_b <- permutest(dbrda(foB, data = xc), permutations = 1000) # pure b

  rm("mems", envir = environment(fit_beta_div))

  # results ----
  return(
    list(
      "component_summary" = ms,
      "relative_contributions" = r_contrib,
      "selected_mems" = mems_string,
      "variation_partitioning" = vp,
      "fraction_abc" = perm_abc,
      "fraction_a" = perm_a,
      "fraction_b" = perm_b,
      "fraction_pure_a" = perm_pure_a,
      "fraction_pure_b" = perm_pure_b
    )
  )
}
