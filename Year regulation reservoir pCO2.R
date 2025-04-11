
rm(list=ls())
library('viridis')


# Define a list of parameter sets for different conditions (wet and dry seasons with or without carbonate buffering)
parameter_sets <- list(
  list(    # Inflow during the wet season
    alk_w = 0.002175,  # Alkalinity in mol/L
    pco2_w = 1957.923769,  # pCO2 in μatm
    kh = 0.038552,  # kh = [CO2(g)]/[CO2(aq)]
    ka1 = 10^-6.376562,  # ka1 = [H]*[HCO3]/[CO2(aq)]
    ka2 = 10^-10.371051  # ka2 = [H]*[CO3]/[HCO3]
  ),
  list(    # Hypolimnion during the wet season
    alk_w = 0.002349,
    pco2_w = 4555.457013,
    kh = 0.047964,
    ka1 = 10^-6.433962,
    ka2 = 10^-10.44611179
  ),
  list(    # Released water during the wet season
    alk_w = 0.002543,
    pco2_w = 3980.487129,
    kh = 0.043993,
    ka1 = 10^-6.410769,
    ka2 = 10^-10.416819 
  ),
  list(    # Inflow without carbonate buffering during the wet season
    alk_w = 0,
    pco2_w = 1957.923769,
    kh = 0.038552,
    ka1 = 10^-6.376562,
    ka2 = 10^-10.371051 
  ),
  list(    # Hypolimnion without carbonate buffering during the wet season
    alk_w = 0,
    pco2_w = 4555.457013,
    kh = 0.047964,
    ka1 = 10^-6.433962,
    ka2 = 10^-10.44611179
  ),
  list(    # Released water without carbonate buffering during the wet season
    alk_w = 0,
    pco2_w = 3980.487129,
    kh = 0.043993,
    ka1 = 10^-6.410769,
    ka2 = 10^-10.416819 
  ),
  
  list(    # Inflow during the dry season
    alk_w = 0.002168,
    pco2_w = 1192.584477,
    kh = 0.044103,
    ka1 = 10^-6.410718,
    ka2 = 10^-10.416170
  ),
  list(    # Hypolimnion during the dry season
    alk_w = 0.002210,
    pco2_w = 3715.870812,
    kh = 0.045579,
    ka1 = 10^-6.419971,
    ka2 = 10^-10.428384
  ),
  list(    # Released water during the dry season
    alk_w = 0.002171,
    pco2_w = 3260.698499,
    kh = 0.045724,
    ka1 = 10^-6.420997,
    ka2 = 10^-10.429814
  ),
  list(    # Inflow without carbonate buffering during the dry season
    alk_w = 0,
    pco2_w = 1192.584477,
    kh = 0.044103,
    ka1 = 10^-6.410718,
    ka2 = 10^-10.416170
  ),
  list(    # Hypolimnion without carbonate buffering during the dry season
    alk_w = 0,
    pco2_w = 3715.870812,
    kh = 0.045579,
    ka1 = 10^-6.419971,
    ka2 = 10^-10.428384
  ),
  list(    # Released water without carbonate buffering during the dry season
    alk_w = 0,
    pco2_w = 3260.698499,
    kh = 0.045724,
    ka1 = 10^-6.420997,
    ka2 = 10^-10.429814
  )
)

## pCO2 in the atmosphere (ppm) 
pco2_atm <- 410 

# Air-water transfer rate constant 
k <- 0.05 

# Maximum simulation time and number of time steps 
t_max <- 100000
t_length <- 100000 #### number of timesteps
t <- seq(0.1, t_max, length.out = t_length)  # Create a time vector

# Function to initialize a matrix of given size, filled with zeros
initialize_matrix <- function(nrow, ncol) {
  matrix(0, nrow = nrow, ncol = ncol)
}

# List to store the simulation results for each parameter set
data_list <- list()

# Loop through each parameter set to run simulations for different conditions
for (i in 1:length(parameter_sets)) {
  # Extract values for each parameter from the current set
  alk_w <- parameter_sets[[i]]$alk_w
  pco2_w <- parameter_sets[[i]]$pco2_w
  kh <- parameter_sets[[i]]$kh
  ka1 <- parameter_sets[[i]]$ka1
  ka2 <- parameter_sets[[i]]$ka2
  
  #### CO2(aq) of water in equilibrium with atmosphere
  co2_atm <- pco2_atm * 10^-6 * kh
  co2_w <- pco2_w * 10^-6 * kh
  
  # Initialize matrices to store the concentrations of CO2, pH, bicarbonate (HCO3-), carbonate (CO32-), DIC, and CO2 degassing rate
  co2 <- initialize_matrix(t_length, 1)
  h <- initialize_matrix(t_length, 1)
  hco3 <- initialize_matrix(t_length, 1)
  co3 <- initialize_matrix(t_length, 1)
  dic <- initialize_matrix(t_length, 1)
  degas_co2 <- initialize_matrix(t_length, 1)
  
  # Define a function to calculate the pH of the system by solving chemical equilibria
  speciate <- function(x, co2_w, alk_w_j, ka1, ka2, kw) {
    abs(-alk_w_j * x^2 + ka1 * co2_w * x + 2 * ka1 * ka2 * co2_w - x^3 + kw * x)
  }
  
  # Define a function to calculate alkalinity based on DIC concentration and pH
  alk_test <- function(x, dic, alk_w_j, ka1, ka2, kw) {
    abs(dic / (1 + x / ka1 + ka2 / x) + 2 * dic / (1 + x / ka2 + x^2 / (ka1 * ka2)) + kw / x - x - alk_w_j)
  }
  
  # Set initial values for pH and other chemical species at the first time step
  h[1, 1] <- optimize(speciate, lower = 10^-9, upper = 10^-4, tol = 10^-18, 
                      co2_w = co2_w, alk_w_j = alk_w, ka1 = ka1, ka2 = ka2, kw = 10^-14)$minimum
  co2[1, 1] <- co2_w
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
                        dic = dic[j, 1], alk_w_j = alk_w, ka1 = ka1, ka2 = ka2, kw = 10^-14)$minimum
    
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
    HCO3 = hco3,   # Bicarbonate concentration
    degas_co2 = degas_co2,  # CO2 degassing rate
    pCO2 = co2 * 10^6 / kh  # pCO2 concentration
  )
}

# Combine all the data frames from the different simulations into one
combined_data <- do.call(cbind, data_list)

# Write the combined data to a CSV file for further analysis or plotting
write.csv(combined_data, file = "Year regulation reservoir.csv", row.names = FALSE)

