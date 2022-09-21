################################################################################
## Functions for setting up the simulation                                    ##
################################################################################

logit <- function(p) {
  log(p/(1-p))
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}


simulate_births <- function(prob.birth, no.mothers) {
  mothers <- matrix(0, nrow=no.mothers, ncol=length(prob.birth))
  # Simulates the birth year
  for(index in 1:length(prob.birth)) {
    mothers[,index] <- rbinom(no.mothers, 1, prob.birth[index])
  }
  return(mothers)
}


create_prob.death.expanded <- function(prob.death, years, periods) {
  # periods should be a vector of length years that maps each year to a period
  # Creates a matrix such that each column is for a specific cohort, 
  #     the probability of dying within the next year
  prob.death.expanded <- matrix(NA, nrow=nrow(prob.death), ncol=nrow(prob.death))
  for(i in 1:length(years)) {
    prob.death.expanded[1:(length(years) -i + 1), i] <-
      sapply(1:(length(years)-i+1), function(x){
        prob.death[x, periods[i+x-1]]
        })
  }
  return(prob.death.expanded)
}

simulate_deaths <- function(birth.year, y.s, prob.death.expanded, 
                            possible.birth.years) {
  # Simulate the death year for one child
  relative.birth.year <- which(possible.birth.years == birth.year)
  death <- F
  i <- 1
  for(year in (birth.year+1):y.s) {
    death <- sample(c(T, F), 1, 
                    prob = c(prob.death.expanded[i, relative.birth.year], 
                             1- prob.death.expanded[i, relative.birth.year]))
    i <- i+1
    if(death == T) {
      death.year <- year
      break
    }
  }
  if(!death) {
    death.year <- NA
  }
  return(death.year)
}

create_sim <- function(m.s, prob.birth, prob.death, y.s, years, periods) {
  # m.s is the age we survey the woman (length is the number of women we survey)
  # prob.death should come in as a matrix of 1q0, 1q1, 1q2, ... for each period
  # prob.birth should be the generic probability of giving birth 
  #     (not conditioned on only having one child) (length is length(15:49) for each age)
  # y.s is the year the survey was conducted
  # This function will remove any women that did not have a birth in the 
  #     observed time span
  
  # returns a data frame with one row per child, containing the following columns:
  # - id: mother's id number
  # - m.s: mother's age at time of survey
  # - m.b: mother's age when they gave birth
  # - y.s: year of survey
  # - y.b: year when they gave birth
  # - death: either NA if child was still alive at time of survey, or year in which they died
  no.women <- length(m.s)
  births <- simulate_births(prob.birth, no.women)
  women.births <- apply(births, 1, function(x){which(x!=0)}) # list of when women had children
  women.validbirths <- lapply(1:length(women.births), function(x) {
    women.births[[x]][which(women.births[[x]] < m.s[x])]
  })
  no.births <- sum(sapply(women.validbirths, length))
  birth.df <- data.frame(id=rep(NA, no.births),
                         m.s=rep(NA, no.births),
                         m.b=rep(NA, no.births))
  i <- 1
  for(woman in 1:no.women){
    if(length(women.validbirths[[woman]])>0) {
      for(birth in 1:length(women.validbirths[[woman]])) {
        birth.df$id[i] <- woman
        birth.df$m.s[i] <- m.s[woman]
        birth.df$m.b[i] <- women.validbirths[[woman]][birth]
        i <- i+1
      }
    }
  }
  birth.df$y.s <- y.s
  birth.df$y.b <- birth.df$y.s - birth.df$m.s + birth.df$m.b
  
  prob.death.expanded <- create_prob.death.expanded(prob.death, years, periods)
  
  birth.df$death <- sapply(1:nrow(birth.df), function(x) {
    simulate_deaths(birth.df$y.b[x], birth.df$y.s[x], 
                    prob.death.expanded=prob.death.expanded, possible.birth.years=years)})
  return(birth.df)
}

# creates a dataframe from full birth history data with the following columns:
# - years: year
# - age: age of children
# - agegp: age group (i.e. [0,1), [1,5), etc.)
# - N: total 
# - Y: deaths 
overallHaz.shorter <- function(fbh) {
  #   year of birth to year of death or year of survey - 1 (i.e., assuming they 
  #   are born/die at the start of the year)
  
  tmp <- lapply(1:nrow(fbh), function(x) {
    create_entries.shorter(fbh[x,])})
  super.df <- data.table::rbindlist(tmp)
  
  suppressMessages(ret <- super.df %>% group_by(years, age, agegp) %>% #change if we want others
                     summarise(N=sum(N), Y=sum(Y)))
  return(ret)
}

create_entries.shorter <- function(fbh.row) {
  max.age <- 49
  y.b <- fbh.row$y.b
  y.s <- fbh.row$y.s
  
  death <- fbh.row$death
  df <- data.frame(years=y.b:(min(y.s,death,na.rm=T)))
  df$age <- 1:nrow(df) - 1
  df$N <- fbh.row$count
  df$Y <- 0
  df$agegp <- cut(df$age, c(0,1,5,max.age), include.lowest = T, right=F)
  df$y.s <- F
  if(fbh.row$y.b != fbh.row$y.s) {
    df <- df[-nrow(df),]
  } else {
    df$y.s <- T
  }
  # note oddness when death > year of survey (this will assume the death occurs at time of survey)
  df$Y[nrow(df)] <- ifelse(!is.na(fbh.row$death),fbh.row$count,0)
  suppressMessages(ret <- data.frame(df %>% group_by(years, age, agegp) %>% 
                                summarise(N=sum(N), Y=sum(Y))))
  return(ret)
}

