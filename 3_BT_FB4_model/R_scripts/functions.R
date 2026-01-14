#functions for BT bioenergetics model

#load packages
library(tidyverse)

#setup grow function 
grow <- function(data, meta, ndays, p0, stress_df,
                 oxycal, EDP) {
  # day-by-day schedules
  p_vec    <- rep_len(p0, ndays)
  epoc_vec <- rep_len(0,  ndays)   
  
  p_vec[stress_df$day]    <- stress_df$p_stress 
  epoc_vec[stress_df$day] <- stress_df$EPOC.J.g
  

  for (i in seq_len(ndays)) {
    
    # ----- Consumption thermal performance (ftC) -----
    data$CG1[i] <- (1/(meta$CTO - meta$CQ)) * log((0.98*(1 - meta$CK1))/(meta$CK1*0.02))
    data$L1[i]  <- exp(data$CG1[i] * (data$temp[i] - meta$CQ))
    data$KA[i]  <- (meta$CK1 * data$L1[i]) / (1 + meta$CK1 * (data$L1[i] - 1))
    
    data$CG2[i] <- (1/(meta$CTL - meta$CTM)) * log((0.98*(1 - meta$CK4))/(meta$CK4*0.02))
    data$L2[i]  <- exp(data$CG2[i] * (meta$CTL - data$temp[i]))
    data$KB[i]  <- (meta$CK4 * data$L2[i]) / (1 + meta$CK4 * (data$L2[i] - 1))
    
    data$ftC[i]  <- data$KA[i] * data$KB[i]
    data$Cmax[i] <- meta$CA * (data$weight[i])^meta$CB
    
    # ----- Use per-day p from schedule -----
    data$Cons.p[i] <- p_vec[i]
    
    data$C[i]      <- data$Cmax[i] * data$Cons.p[i] * data$ftC[i]     # proportion of Cmax * temp effect
    data$Cons.g[i] <- data$C[i] * data$weight[i]                       
    data$Cons.J[i] <- data$Cons.g[i] * EDP                              # J ingested
    
    # ----- Wastes -----
    data$Eg[i] <- meta$FA * (data$temp[i])^(meta$FB) * exp(meta$FG * data$Cons.p[i]) * data$Cons.J[i]
    data$Ex[i] <- meta$UA * (data$temp[i])^(meta$UB) * exp(meta$UG * data$Cons.p[i]) * (data$Cons.J[i] - data$Eg[i])
    
    # ----- SDA -----
    data$SDA[i] <- meta$SDA * (data$Cons.J[i] - data$Eg[i])
    
    # ----- Metabolism thermal performance (ftR) -----
    data$Z[i] <- log(meta$RQ) * (meta$RTM - meta$RTO)
    data$Y[i] <- log(meta$RQ) * (meta$RTM - meta$RTO + 2)
    data$X[i] <- (data$Z[i]^2 * (1 + sqrt(1 + 40/data$Y[i]))^2) / 400
    data$V[i] <- (meta$RTM - data$temp[i]) / (meta$RTM - meta$RTO)
    
    data$ftR[i]  <- ifelse(data$temp[i] < meta$RTM,
                           data$V[i]^data$X[i] * exp(data$X[i] * (1 - data$V[i])),
                           1e-6)
    data$Rmax[i] <- meta$RA * (data$weight[i])^meta$RB
    data$Met.J[i] <- data$Rmax[i] * data$ftR[i] * data$weight[i] * oxycal
    
    # ----- EPOC (one-time cost today) -----
    data$EPOC.J[i] <- epoc_vec[i] * data$weight[i]
    
    # ----- Growth  -----
    data$Growth.J[i] <- data$Cons.J[i] - (data$Met.J[i] + data$SDA[i] + data$Eg[i] + data$Ex[i] + data$EPOC.J[i])
    data$Growth.g[i] <- data$Growth.J[i] / meta$ED
    
    # ----- update weight -----
      data$E[i + 1]      <- data$E[i] + data$Growth.J[i]
      data$weight[i + 1] <- data$E[i + 1] / meta$ED
    }
  
  list(data = dplyr::filter(data, day <= ndays),
       final_weight = data$weight[ndays])
}

#################################################################
#Fit-p fxn
#Determine p0 based on initial and final weights using bioenergetics model
#Uses bisection method to find p0 that achieves target final weight over 120 days

fit.p <- function(initial_g, final_g, sim_template, meta, oxycal, EDP,
                  W.tol = 0.5, max.iter = 100) {
  # initial_g: Initial weight in grams
  # final_g: Target final weight in grams
  # sim_template: Dataframe with temperature profile (day, temp columns)
  # meta: Brook trout bioenergetics parameters (bt)
  # oxycal: Oxygen calorie conversion (13560)
  # EDP: Prey energy density (4000)
  # W.tol: Weight tolerance for convergence (default 0.5g)
  # max.iter: Maximum iterations (default 100)

  # Initialize bisection search
  n.iter <- 0
  p.max  <- 1
  p.min  <- 0
  p      <- 0.5  # Initial guess

  # Create empty stress_df (no stress events for baseline growth)
  stress_df <- data.frame(
    day = integer(0),
    p_stress = numeric(0),
    EPOC.J.g = numeric(0)
  )

  # Create simulation dataframe (copy template to avoid modification)
  sim <- sim_template
  sim$weight <- NA
  sim$E <- NA
  sim$weight[1] <- initial_g
  sim$E[1] <- initial_g * meta$ED

  # Initial run
  output <- grow(data = sim, meta = meta, ndays = 120, p0 = p,
                 stress_df = stress_df, oxycal = oxycal, EDP = EDP)
  W.p <- output[[2]]

  # Bisection loop to find p0 that produces target final weight
  while((n.iter <= max.iter) & (abs(W.p - final_g) > W.tol)) {
    n.iter <- n.iter + 1

    # Adjust search bounds
    if(W.p > final_g) {
      p.max <- p
    } else {
      p.min <- p
    }

    # Update p to midpoint
    p <- (p.min + p.max) / 2

    # Reset simulation
    sim <- sim_template
    sim$weight <- NA
    sim$E <- NA
    sim$weight[1] <- initial_g
    sim$E[1] <- initial_g * meta$ED

    # Run grow with new p
    output <- grow(data = sim, meta = meta, ndays = 120, p0 = p,
                   stress_df = stress_df, oxycal = oxycal, EDP = EDP)
    W.p <- output[[2]]
  }

  return(p)
}

#################################################################
#setup fxn for 72 scenarios based on (1) number of C&R events, (2) temp, (3) weight, (4) EPOC level

run_scenario <- function(scenario_row, sim_template, meta, ndays, p0, y_spacing, oxycal, EDP) {
  # scenario_row: single row from stress parameter grid with columns:
  #   Temperature..C, Weight_Class, Weight_kg, F_after, EPOC_level, EPOC.J.g, n_events

  # Set up simulation data frame (copy template)
  sim <- sim_template

  # Set initial weight from Weight_kg (convert kg to g)
  sim$weight[1] <- scenario_row$Weight_kg * 1000
  sim$E[1] <- sim$weight[1] * meta$ED

  # Build stress_df for this scenario
  if (scenario_row$n_events > 0) {
    stress_df <- scenario_row %>%
      slice(rep(1, scenario_row$n_events)) %>%
      mutate(
        day = row_number() * y_spacing,
        p_stress = F_after * p0
      )
  } else {
    # No events: create empty dataframe with required columns
    stress_df <- data.frame(
      day = integer(0),
      p_stress = numeric(0),
      EPOC.J.g = numeric(0)
    )
  }

  # Run grow function
  output <- grow(
    data = sim,
    meta = meta,
    ndays = ndays,
    p0 = p0,
    stress_df = stress_df,
    oxycal = oxycal,
    EDP = EDP
  )

  # Extract final weight and energetics
  energetics <- output[[1]]
  final_weight <- output[[2]]
  
  event_day <- y_spacing
  
  # Extract Net_Energy (Growth.J) for 0 and 1 event scenarios
  if (scenario_row$n_events == 0) {
    net_energy <- energetics$Growth.J[energetics$day == event_day]
  } else if (scenario_row$n_events == 1) {
    net_energy <- energetics$Growth.J[energetics$day == y_spacing]
  } else {
    net_energy <- NA
  }

  # Return summary
  data.frame(
    EPOC_level = scenario_row$EPOC_level,
    Temperature = scenario_row$Temperature..C,
    Weight_Class = scenario_row$Weight_Class,
    n_events = scenario_row$n_events,
    initial_weight_g = sim$weight[1],
    final_weight_g = final_weight,
    growth_g = final_weight - sim$weight[1],
    percent_growth = ((final_weight - sim$weight[1]) / sim$weight[1]) * 100,
    Net_Energy = net_energy
  )
}

. 


