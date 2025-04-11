
rm(list = ls())
library('viridis')

# Define a list of parameter sets for different conditions (wet and dry seasons with or without carbonate buffering)
parameter_sets <- list(
  list(  # Inflow during the wet season
    alk_gw = 0.002375,  # Alkalinity in mol/L
    pco2_gw = 1424.458689, # pCO2 in μatm
    kh = 0.037757,  # kh = [CO2(g)]/[CO2(aq)]
    ka1 = 10^-6.371864,  # ka1 = [H]*[HCO3]/[CO2(aq)]
    ka2 = 10^-10.365390    # ka2 = [H]*[CO3]/[HCO3]
  ),
  list( # Hypolimnion during the wet season
    alk_gw = 0.002470,
    pco2_gw = 4521.272808,
    kh = 0.043237,
    ka1 = 10^-6.405798,
    ka2 = 10^-10.410097 
  ),
  list( # Released water during the wet season
    alk_gw = 0.002401  ,
    pco2_gw = 3352.313510  ,
    kh = 0.040469  ,
    ka1 = 10^-6.389039  ,
    ka2 = 10^-10.388382 
  ),
  list( # Inflow without carbonate buffering during the wet season
    alk_gw = 0,
    pco2_gw = 1424.458689  ,
    kh = 0.037757  ,
    ka1 = 10^-6.371864  ,
    ka2 = 10^-10.365390 
  ),
  list( # Hypolimnion without carbonate buffering during the wet season
    alk_gw = 0,
    pco2_gw = 4521.272808  ,
    kh = 0.043237 ,
    ka1 = 10^-6.405798  ,
    ka2 = 10^-10.410097 
  ),
  list( # Released water without carbonate buffering during the wet season
    alk_gw = 0 ,
    pco2_gw = 3352.313510  ,
    kh = 0.040469  ,
    ka1 = 10^-6.389039  ,
    ka2 = 10^-10.388382 
  ),
  list( # Inflow during the dry season
    alk_gw = 0.002702  ,
    pco2_gw = 1285.382132  ,
    kh =0.044925  ,
    ka1 = 10^-6.416424  ,
    ka2 = 10^-10.424114 
  ),
  list(  # Hypolimnion during the dry season
    alk_gw = 0.002416  ,
    pco2_gw = 3357.912385  ,
    kh = 0.042705  ,
    ka1 = 10^-6.402887  ,
    ka2 = 10^-10.406573 
  ),
  list( # Released water during the dry season
    alk_gw = 0.002365  ,
    pco2_gw = 2812.471512  ,
    kh = 0.042408  ,
    ka1 = 10^-6.401182  ,
    ka2 = 10^-10.404460 
  ),
  list( # Inflow without carbonate buffering during the dry season
    alk_gw =0,
    pco2_gw = 1285.382132  ,
    kh =0.044925  ,
    ka1 = 10^-6.416424  ,
    ka2 = 10^-10.424114 
  ),
  list( # Hypolimnion without carbonate buffering during the dry season
    alk_gw = 0,
    pco2_gw = 3357.912385  ,
    kh = 0.042705  ,
    ka1 = 10^-6.402887  ,
    ka2 = 10^-10.406573 
  ),
  list( # Released water without carbonate buffering during the dry season
    alk_gw = 0 ,
    pco2_gw = 2812.471512  ,
    kh = 0.042408  ,
    ka1 = 10^-6.401182  ,
    ka2 = 10^-10.404460 
  )
)


pco2_atm <- 410  # pCO2 in the atmosphere (ppm)  
k <- 0.05    # Air-water transfer rate constant 
t_max <- 100000   # Maximum simulation time and number of time steps
t_length <- 100000   # number of timesteps
t <- seq(0.1, t_max, length.out = t_length)   # Create a time vector

# Function to initialize a matrix of given size, filled with zeros
initialize_matrix <- function(nrow, ncol) {
  matrix(0, nrow = nrow, ncol = ncol)
}

# List to store the simulation results for each parameter set
data_list <- list()

# Loop through each parameter set to run simulations for different conditions
for (i in 1:length(parameter_sets)) {
  # Extract values for each parameter from the current set
  alk_gw <- parameter_sets[[i]]$alk_gw
  pco2_gw <- parameter_sets[[i]]$pco2_gw
  kh <- parameter_sets[[i]]$kh
  ka1 <- parameter_sets[[i]]$ka1
  ka2 <- parameter_sets[[i]]$ka2
 
   #### CO2(aq) of water in equilibrium with atmosphere
  co2_atm <- pco2_atm * 10^-6 * kh
  co2_gw <- pco2_gw * 10^-6 * kh
  
  # Initialize matrices to store the concentrations of CO2, pH, bicarbonate (HCO3-), carbonate (CO32-), DIC, and CO2 degassing rate
  co2 <- initialize_matrix(t_length, 1)
  h <- initialize_matrix(t_length, 1)
  hco3 <- initialize_matrix(t_length, 1)
  co3 <- initialize_matrix(t_length, 1)
  dic <- initialize_matrix(t_length, 1)
  degas_co2 <- initialize_matrix(t_length, 1)
  
  # Define a function to calculate the pH of the system by solving chemical equilibria
  speciate <- function(x, co2_gw, alk_gw_j, ka1, ka2, kw) {
    abs(-alk_gw_j * x^2 + ka1 * co2_gw * x + 2 * ka1 * ka2 * co2_gw - x^3 + kw * x)
  }
  
  # Define a function to calculate alkalinity based on DIC concentration and pH
  alk_test <- function(x, dic, alk_gw_j, ka1, ka2, kw) {
    abs(dic / (1 + x / ka1 + ka2 / x) + 2 * dic / (1 + x / ka2 + x^2 / (ka1 * ka2)) + kw / x - x - alk_gw_j)
  }
  
  # Set initial values for pH and other chemical species at the first time step
  h[1, 1] <- optimize(speciate, lower = 10^-9, upper = 10^-4, tol = 10^-18, 
                      co2_gw = co2_gw, alk_gw_j = alk_gw, ka1 = ka1, ka2 = ka2, kw = 10^-14)$minimum
  co2[1, 1] <- co2_gw
  hco3[1, 1] <- ka1 * co2[1, 1] / h[1, 1]
  co3[1, 1] <- ka2 * hco3[1, 1] / h[1, 1]
  dic[1, 1] <- co2[1, 1] + hco3[1, 1] + co3[1, 1]
  
  # Loop through each time step to calculate changes in pH and species concentrations
  for (j in 2:t_length) {
    # Calculate the rate of CO2 degassing based on the difference between water CO2 and atmospheric CO2
    degas_co2[j, 1] <- (t_max / t_length) * k * (co2[j - 1, 1] - co2_atm)
    
    # Update the DIC concentration by subtracting the degassed CO2 amount
    dic[j, 1] <- dic[j - 1, 1] - degas_co2[j, 1]
    
    # Calculate pH by optimizing the alkalinity equation at the current DIC concentration
    h[j, 1] <- optimize(alk_test, lower = 10^-12, upper = 10^-3, tol = 10^-14, 
                        dic = dic[j, 1], alk_gw_j = alk_gw, ka1 = ka1, ka2 = ka2, kw = 10^-14)$minimum
    
    # Calculate concentrations of CO2, HCO3-, and CO32- based on pH
    co2[j, 1] <- dic[j, 1] / (1 + ka1 / h[j, 1] + ka1 * ka2 / h[j, 1]^2)
    hco3[j, 1] <- ka1 * co2[j, 1] / h[j, 1]
    co3[j, 1] <- ka2 * hco3[j, 1] / h[j, 1]
  }
  
  # Store the simulation results for the current parameter set in the data list
  data_list[[i]] <- data.frame(
    Time = t * k,  # Time, scaled by the transfer rate constant
    pH = -log10(h),  # pH, calculated as the negative logarithm of H+ concentration
    DIC = dic,  # Dissolved Inorganic Carbon concentration
    CO2 = co2,  # CO2 concentration
    CO3 = co3,  # Carbonate concentration
    HCO3 = hco3,  # Bicarbonate concentration
    degas_co2 = degas_co2,  # CO2 degassing rate
    pCO2 = co2 * 10^6 / kh  # pCO2 concentration
  )
}

# Combine all the data frames from the different simulations into one
combined_data <- do.call(cbind, data_list)

# Write the combined data to a CSV file for further analysis or plotting
write.csv(combined_data, file = "Seasonal regulation reservoir.csv", row.names = FALSE)


