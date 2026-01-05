# ================================================================
# Classical Analysis of Crime Data (RA 9165 Cases, 2013–2024)
# Using Descriptive and Poisson-Based Time Series Models
# Author: [Your Name]
# ================================================================

# ---- Load Required Libraries ----
library(ggplot2)       # For visualization
library(dplyr)         # For data transformation
library(lubridate)     # For parsing date formats
library(tscount)       # For time series GLMs with count data
library(forecast)      # For residual diagnostics

# ---- Step 1: Load and Prepare Dataset ----
# Load monthly drug-related case counts under RA 9165 from 2013–2024

data <- read.csv("C:/Users/ACER/Desktop/DC.csv")
colnames(data) <- c("Date", "Cases")
data$Date <- dmy_hm(data$Date)  # Convert to date format
data <- data %>% arrange(Date)

# ---- Step 2: Traditional Descriptive Analysis ----
# Law enforcement typically relies on summaries and visual patterns

# Summary statistics
summary_stats <- summary(data$Cases)
sd_cases <- sd(data$Cases)

cat("Summary of Monthly Drug Cases (RA 9165):\n")
print(summary_stats)
cat("Standard Deviation:", round(sd_cases, 2), "\n")

# Add Index column
data$Index <- 1:nrow(data)

# Identify top 3 highest points
top3 <- data %>% 
  arrange(desc(Cases)) %>% 
  slice(1:3)

# Flag column for top 3
data$IsTop3 <- data$Index %in% top3$Index

# Plot
ggplot(data, aes(x = Index, y = Cases)) +
  geom_line(color = "royalblue", linewidth = 1.5) +  # thicker & more vivid
  geom_point(color = "darkred") +
  geom_text(aes(label = Index), 
            vjust = -0.7, 
            size = ifelse(data$IsTop3, 5, 2.5), 
            color = "black") +
  labs(title = "Monthly Reported Drug Cases (RA 9165): Jan 2013 – Dec 2024",
       x = "Observation Index", y = "Number of Cases") +
  theme_minimal() +
  scale_x_continuous(breaks = data$Index) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# Histogram of Count Distribution
ggplot(data, aes(x = Cases)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Monthly Drug Cases",
       x = "Monthly Case Count", y = "Frequency") +
  theme_minimal()

# ---- Step 3: Poisson Autoregressive (PAR(1)) Model ----
# Models the current count as dependent on the previous count

par_model <- tsglm(data$Cases, model = list(past_obs = 1), distr = "poisson")

cat("\n--- Poisson AR(1) Model Summary ---\n")
summary(par_model)

# ---- Step 4: Poisson-ARCH Model (AR + Lagged Mean) ----
# Adds lagged conditional mean to account for volatility in counts

par_arch_model <- tsglm(data$Cases, model = list(past_obs = 1, past_mean = 1), distr = "poisson")

cat("\n--- Poisson-ARCH(1,1) Model Summary ---\n")
summary(par_arch_model)

# ---- Step 5: Poisson GLM with Lagged Predictor ----
# Simple GLM assuming log link and one lagged observation

data$Lag1 <- lag(data$Cases, 1)
glm_data <- na.omit(data)

glm_model <- glm(Cases ~ Lag1, data = glm_data, family = "poisson")

cat("\n--- Poisson GLM Summary ---\n")
summary(glm_model)

# ---- Step 6: Check for Overdispersion in GLM ----
# Compare residual deviance to degrees of freedom

cat("\n--- Overdispersion Check ---\n")
cat("Residual Deviance:", deviance(glm_model), "\n")
cat("Degrees of Freedom:", df.residual(glm_model), "\n")

# ---- Step 7: Interpretation and Transition Notes ----
# Observations:
# - Time plot reveals spikes in early 2020 and late 2023.
# - Histogram suggests right-skew and overdispersion.
# - Poisson AR captures short-term memory but not structural breaks.
# - Poisson-ARCH improves AIC slightly but remains limited.
# - GLM shows residual deviance >> degrees of freedom ⇒ overdispersion.
#
# Conclusion:
# Classical Poisson-based models provide baseline insight but fail to
# accommodate overdispersion, time-varying effects, or interventions.
#
# Recommendation:
# Transition to Bayesian Integer-Valued Transfer Function Models (TFMs)
# using Generalized Poisson or Negative Binomial distributions.
# These allow for:
# - Modeling intervention timing and duration
# - Posterior uncertainty quantification
# - Model comparison via DIC instead of AIC
# - Proper treatment of overdispersed, real-world count data

# ---- Next Step: Bayesian TFM Modeling ----
# Proceed to:
# - GP-TFM1, GP-TFM2 (Generalized Poisson)
# - NB-TFM1, NB-TFM2 (Negative Binomial)
# using MCMC simulation with prior/posterior diagnostics
