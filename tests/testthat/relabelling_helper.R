library(partitions)

paste_pattern <- function(x) paste0(x, collapse="")

match_partitions <- function(patterns, ref_patterns){
  temp <- lapply(strsplit(as.character(patterns), ""), as.numeric)
  res <- rep(NA, length(temp))
  for(i in 1:length(temp)){
    if(paste_pattern(temp[[i]]) %in% ref_patterns){
      res[i] <- paste_pattern(temp[[i]])
    }else{
      perm <- t(perms(max(temp[[i]])))
      for(row in 1:nrow(perm)){
        if(paste_pattern(perm[row, ][temp[[i]]]) %in% ref_patterns){
          # cat(paste_pattern(temp[[i]]), " -- ", perm[row, ], " -- ", paste_pattern(perm[row, ][temp[[i]]]), "\n")
          res[i] <- paste_pattern(perm[row, ][temp[[i]]])
        }
      }
    }
  }
  res
}

helper_summarise_partitions <- function(clustering){
  patterns_obs <- apply(clustering, 1, paste_pattern)
  df1 <- data.frame(pattern = patterns_obs) %>%
    group_by(pattern) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(freq = n / sum(n)) %>%
    mutate(pattern = as.character(pattern))


  df_empirical <- df1 %>%
    mutate(matched_patterns = match_partitions(df1$pattern, df_correct$pattern)) %>%
    group_by(matched_patterns) %>%
    summarise(count = sum(n)) %>%
    ungroup() %>%
    mutate(freq = count / sum(count)) %>%
    rename(pattern = matched_patterns)

  df_empirical
}
