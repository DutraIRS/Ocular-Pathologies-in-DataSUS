prepare_data <- function(disease="glaucoma") {
  # Read data
  data <- read.csv("data/sus_ocular_data.csv")
  
  # -1 is na
  data[data == -1] <- NA
  
  # Convert categorical variables to factors
  data$state <- as.factor(data$state)
  
  # Age in relation to 60 (effect in decades because poisson was exploding with e^85)
  data$age_relative <- (data$age - 60) / 10
  
  # Remove 2025
  data <- subset(data, year != 2025)
  
  # Year in relation to start date
  data$year_relative <- data$year - min(data$year)
  
  # Use log pop
  data$log_pop <- log(data$population)
  
  # Generic disease
  data$disease <- data[[disease]]
  
  data <- data[c("state", "year_relative", "age_relative", "log_pop", "disease")]
  data <- na.omit(data)
  
  data$state_id <- as.integer(as.factor(data$state))
  data$year_id <- data$year_relative + 1
  
  return(data)
}