
rm(list=ls())
library('viridis')


## Set the pCO2 and delta-13C values for atmosphere
pco2_atm <- 410 ## pCO2 in the atmosphere (ppm)
delta_atm <- -9  # Atmospheric delta-13C value
## PDB ratio for carbon isotopes (used to convert delta-13C)
r_pdb <- 0.011179

# Define a list of parameter sets for different conditions (wet and dry seasons with or without carbonate buffering)
params <- list(  #Inflow during the wet season
  list(alk_gw = 0.002375,  # Alkalinity in mol/L
       pco2_gw =1424.458689,   # pCO2 in ppm
       delta_dic_o = -9.295000, # Water delta-13C
       alpha = 1,   # kinetic fractionation factor during degassing
       kh =0.037757,  # kh = [CO2(g)]/[CO2(aq)]
       ka1 <- 10^-6.371864,
       ka1_13 <- 10^-6.352424, # ka1 = [H]*[HCO3]/[CO2(aq)]
       ka2 <- 10^-10.365390 ,
       ka2_13 <- 10^-10.353760),  # ka2 = [H]*[CO3]/[HCO3]
  list(alk_gw =0.002470,  #Hypolimnion during the wet season
       pco2_gw =4521.272808,
       delta_dic_o = -9.727933, 
       alpha = 1, 
       kh = 0.043237, 
       ka1 <- 10^-6.405798,
       ka1_13 <- 10^-6.390200 ,
       ka2 <- 10^-10.410097 ,
       ka2_13 <- 10^-10.400781),
  list(alk_gw =0.002401,    #Released water during the wet season
       pco2_gw = 3352.313510, 
       delta_dic_o = 	-9.685750, 
       alpha = 1, 
       kh = 0.040469, 
       ka1 <- 10^-6.389039,
       ka1_13 <- 10^-6.371621,
       ka2 <- 10^-10.388382,
       ka2_13 <- 10^-10.377971),
  list(alk_gw = 0.002702,  #Inflow during the dry season
       pco2_gw =1285.382132, 
       delta_dic_o = -8.684000, 
       alpha = 1, 
       kh =0.044925, 
       ka1 <- 10^-6.416424,
       ka1_13 <- 10^-6.402033,
       ka2 <- 10^-10.424114,
       ka2_13 <- 10^-10.415525 ),
  list(alk_gw =0.002416,   #Hypolimnion during the dry season
       pco2_gw = 3357.912385 , 
       delta_dic_o = -9.908667, 
       alpha = 1, 
       kh = 0.042705, 
       ka1 <- 10^-6.402887,
       ka1_13 <- 10^-6.387025,
       ka2 <- 10^-10.406573,
       ka2_13 <- 10^-10.397098),
  list(alk_gw = 0.002365,   #Released water during the dry season
       pco2_gw =2812.471512 , 
       delta_dic_o = -10.097500, 
       alpha = 1, 
       kh =0.042408, 
       ka1 <- 10^-6.401182,
       ka1_13 <- 10^-6.385155,
       ka2 <- 10^-10.404460,
       ka2_13 <- 10^-10.394887),
  list(alk_gw = 0,  #Inflow without carbonate buffering during the wet season
       pco2_gw =1424.458689, 
       delta_dic_o = -9.295000, 
       alpha = 1, 
       kh =0.037757, 
       ka1 <- 10^-6.371864,
       ka1_13 <- 10^-6.352424,
       ka2 <- 10^-10.365390 ,
       ka2_13 <- 10^-10.353760),
  list(alk_gw =0,   #Hypolimnion without carbonate buffering during the wet season
       pco2_gw =4521.272808, 
       delta_dic_o = -9.727933, 
       alpha = 1, 
       kh = 0.043237, 
       ka1 <- 10^-6.405798,
       ka1_13 <- 10^-6.390200 ,
       ka2 <- 10^-10.410097 ,
       ka2_13 <- 10^-10.400781),
  list(alk_gw =0, #Released water without carbonate buffering during the wet season
       pco2_gw = 3352.313510, 
       delta_dic_o = 	-9.685750, 
       alpha = 1, 
       kh = 0.040469, 
       ka1 <- 10^-6.389039,
       ka1_13 <- 10^-6.371621,
       ka2 <- 10^-10.388382,
       ka2_13 <- 10^-10.377971 ),
  list(alk_gw = 0, #Inflow without carbonate buffering during the dry season
       pco2_gw =1285.382132, 
       delta_dic_o = -8.684000, 
       alpha = 1, 
       kh =0.044925, 
       ka1 <- 10^-6.416424,
       ka1_13 <- 10^-6.402033,
       ka2 <- 10^-10.424114,
       ka2_13 <- 10^-10.415525),
  list(alk_gw =0,  #Hypolimnion without carbonate buffering during the dry season
       pco2_gw = 3357.912385 , 
       delta_dic_o = -9.908667, 
       alpha = 1, 
       kh = 0.042705, 
       ka1 <- 10^-6.402887,
       ka1_13 <- 10^-6.387025,
       ka2 <- 10^-10.406573,
       ka2_13 <- 10^-10.397098),
  list(alk_gw = 0, #Released water without carbonate buffering during the dry season
       pco2_gw =2812.471512 , 
       delta_dic_o = -10.097500, 
       alpha = 1, 
       kh =0.042408, 
       ka1 <- 10^-6.401182,
       ka1_13 <- 10^-6.385155,
       ka2 <- 10^-10.404460,
       ka2_13 <- 10^-10.394887)
)


#### degassing rate
k <- 0.05
#### total simulation time
t_max <- 100000
#### number of timesteps
t_length <- 100000
#### array of timesteps
t <- seq(0.1, t_max, length.out = t_length)

# List to store the simulation results for each parameter set
delta_dic_list <- list()

# Loop through each parameter set to run simulations for different conditions
for (i in seq_along(params)) {
  param_set <- params[[i]]
  
  
  # Extract values for each parameter from the current set 
  alk_gw <- param_set$alk_gw
  pco2_gw <- param_set$pco2_gw
  delta_dic_o <- param_set$delta_dic_o
  alpha <- param_set$alpha
  kh <- param_set$kh
  col <- param_set$col
  
 
  kw <- 10^-14
  
  # Calculate CO2 in equilibrium with the atmosphere and initial carbon isotope values
  co2_atm <- pco2_atm * 10^-6 * kh
  co2_gw <- pco2_gw * 10^-6 * kh
  
  c13_atm <- ((delta_atm / 1000 + 1) * r_pdb) * co2_atm
  
  # Initialize arrays for storing the results at each timestep
  dic <- array(dim = c(length(t), 1))
  co2 <- dic
  hco3 <- dic
  co3 <- dic
  h <- dic
  c13_dic <- dic
  c13_co2 <- dic
  c13_hco3 <- dic
  c13_co3 <- dic
  delta_dic <- dic
  delta_co2 <- dic
  delta_hco3 <- dic
  delta_co3 <- dic
  degas_co2 <- dic
  
  # Initial condition for pH (h) at the first timestep
  h[1] <- optimize(function(x) {abs(-alk_gw * x^2 + ka1 * co2_gw * x + 2 * ka1 * ka2 * co2_gw - x^3 + kw * x)}, lower = 10^-9, upper = 10^-4, tol = 10^-18)$minimum
  co2[1] <- co2_gw
  hco3[1] <- ka1 * co2[1] / h[1]
  co3[1] <- ka2 * hco3[1] / h[1]
  dic[1] <- co2[1] + hco3[1] + co3[1]
  c13_dic[1] <- ((delta_dic_o / 1000 + 1) * r_pdb) * dic[1]
  c13_co2[1] <- c13_dic[1] / (1 + ka1_13 / h[1] + ka1_13 * ka2_13 / h[1]^2)
  c13_hco3[1] <- ka1_13 * c13_co2[1] / h[1]
  c13_co3[1] <- ka2_13 * c13_hco3[1] / h[1]
  
  # Calculate initial delta values for different species
  delta_dic[1] <- delta_dic_o
  delta_co2[1] <- delta_dic[1] - 9 * hco3[1] / dic[1] - 8 * co3[1] / dic[1]
  delta_hco3[1] <- delta_co2[1] + 9
  delta_co3[1] <- delta_co2[1] + 8
  
  
  # Loop through each timestep to simulate the process
  for (j in 2:length(t)) {  
    degas_co2[j] <- (t_max / t_length) * k * (co2[j-1] - co2_atm)    # Degassing rate for CO2
    dic[j] <- dic[j-1] - degas_co2[j]    # Update DIC concentration
    h[j] <- optimize(function(x) {abs(dic[j] / (1 + x / ka1 + ka2 / x) + 2 * dic[j] / (1 + x / ka2 + x^2 / (ka1 * ka2)) + kw / x - x - alk_gw)}, lower = 10^-12, upper = 10^-3, tol = 10^-14)$minimum
    co2[j] <- dic[j] / (1 + ka1 / h[j] + ka1 * ka2 / h[j]^2)
    hco3[j] <- ka1 * co2[j] / h[j]
    co3[j] <- ka2 * hco3[j] / h[j]
    
    # Calculate 13C isotopes for different species
    c_13_stream <- ((delta_co2[j-1] / 1000 + 1) * r_pdb) * co2[j-1]
    c_13_degas <- -alpha * (t_max / t_length) * k * (c_13_stream - c13_atm)
    c_13_dic <- ((delta_dic[j-1] / 1000 + 1) * r_pdb) * dic[j-1]
    c_13_new_dic <- c_13_dic + c_13_degas
    delta_dic[j] <- (c_13_new_dic / dic[j] / r_pdb - 1) * 1000
    delta_co2[j] <- delta_dic[j] - 9 * hco3[j] / dic[j] - 8 * co3[j] / dic[j]
    delta_hco3[j] <- delta_co2[j] + 9
    delta_co3[j] <- delta_co2[j] + 8
    
    c13_dic[j] <- c13_dic[j-1] + c_13_degas
    c13_co2[j] <- c13_dic[j] / (1 + ka1_13 / h[j] + ka1_13 * ka2_13 / h[j]^2)
    c13_hco3[j] <- ka1_13 * c13_co2[j] / h[j]
    c13_co3[j] <- ka2_13 * c13_hco3[j] / h[j]
  }
  
  # Store the delta-13C results for the current parameter set
  delta_dic_list[[i]] <- delta_dic
}

# Combine all the delta-13C data into one matrix for export
delta_dic_all <- do.call(cbind, delta_dic_list)

# Create a data frame for time and results
time_df <- data.frame(Time = t * k)

# Combine time data with delta-13C results for export
combined_data <- cbind(time_df, delta_dic_all)

# Save the combined data to a CSV file
write.csv(combined_data, "Seasonal regulation reservoir 13C.csv", row.names = FALSE)


