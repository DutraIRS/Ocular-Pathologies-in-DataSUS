
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)
library(dplyr)

age_codes <- paste(c(
  93088,  # 25 a 29
  93089,  # 30 a 34
  93090,  # 35 a 39
  93091,  # 40 a 44
  93092,  # 45 a 49
  93093,  # 50 a 54
  93094,  # 55 a 59
  93095,  # 60 a 64
  93096,  # 65 a 69
  93097,  # 70 a 74
  93098,  # 75 a 79
  49108,  # 80 a 84
  49109,  # 85 a 89
  49110   # 90 ou mais
), collapse = ",")

# Year codes: 2010-2060
# 2010=4336, 2011=12037, 2012=13242, then 2013+k -> 49029+k
year_codes <- paste(c(
  4336, 12037, 13242,       # 2010, 2011, 2012
  49029:49076               # 2013 through 2060
), collapse = ",")

api_url <- paste0(
  "https://apisidra.ibge.gov.br/values",
  "/t/7358",          # table
  "/n3/all",          # all states
  "/v/606",           # variable: population
  "/p/all",           # period (reference)
  "/c2/6794",         # sex = total
  "/c287/", age_codes,
  "/c1933/", year_codes
)

cat("Fetching from SIDRA API...\n")
cat("URL length:", nchar(api_url), "chars\n")

response <- tryCatch(
  fromJSON(api_url, flatten = TRUE),
  error = function(e) {
    cat("Direct fetch failed:", conditionMessage(e), "\n")
    cat("Trying in batches...\n")
    NULL
  }
)

if (is.null(response)) {
  # Split year codes into two batches
  batch1 <- paste(c(4336, 12037, 13242, 49029:49052), collapse = ",")  # 2010-2035
  batch2 <- paste(49053:49076, collapse = ",")                         # 2036-2060

  fetch_batch <- function(yr_codes) {
    url <- paste0(
      "https://apisidra.ibge.gov.br/values",
      "/t/7358/n3/all/v/606/p/all",
      "/c2/6794/c287/", age_codes,
      "/c1933/", yr_codes
    )
    fromJSON(url, flatten = TRUE)
  }

  cat("  Batch 1 (2010-2035)...\n")
  r1 <- fetch_batch(batch1)
  cat("  Batch 2 (2036-2060)...\n")
  r2 <- fetch_batch(batch2)

  # Combine (skip header row in batch 2)
  response <- rbind(r1, r2[-1, ])
}

cat("Received", nrow(response) - 1, "data rows.\n")
headers_raw <- as.character(unlist(response[1, ]))
df <- as.data.frame(response[-1, ], stringsAsFactors = FALSE)
names(df) <- make.unique(headers_raw)   # avoids duplicate-name errors

cat("Columns (", ncol(df), "):", paste(names(df), collapse = " | "), "\n")

# Use positional indices (robust to locale-dependent names)
COL_VALUE      <- 5
COL_STATE_CODE <- 6
COL_AGE_LABEL  <- 15
COL_YEAR_LABEL <- 17    # projection year (c1933), NOT the period column

cat("  State code col [", COL_STATE_CODE, "]:", names(df)[COL_STATE_CODE], "\n")
cat("  Age label  col [", COL_AGE_LABEL, "]:", names(df)[COL_AGE_LABEL], "\n")
cat("  Year label col [", COL_YEAR_LABEL, "]:", names(df)[COL_YEAR_LABEL], "\n")
cat("  Value      col [", COL_VALUE, "]:", names(df)[COL_VALUE], "\n")

state_map <- c(
  "11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR", "15" = "PA",
  "16" = "AP", "17" = "TO",
  "21" = "MA", "22" = "PI", "23" = "CE", "24" = "RN", "25" = "PB",
  "26" = "PE", "27" = "AL", "28" = "SE", "29" = "BA",
  "31" = "MG", "32" = "ES", "33" = "RJ", "35" = "SP",
  "41" = "PR", "42" = "SC", "43" = "RS",
  "50" = "MS", "51" = "MT", "52" = "GO", "53" = "DF"
)

# Map 5-year age groups to the project's bins
age_5yr_to_bin <- c(
  "25" = 25, "30" = 25,
  "35" = 35, "40" = 35,
  "45" = 45, "50" = 45,
  "55" = 55, "60" = 55,
  "65" = 65, "70" = 65,
  "75" = 75, "80" = 75, "85" = 75,
  "90" = 90
)

df$age_lower <- as.integer(gsub("^(\\d+).*", "\\1", df[[COL_AGE_LABEL]]))
df$state     <- state_map[as.character(df[[COL_STATE_CODE]])]
df$year      <- as.integer(df[[COL_YEAR_LABEL]])
df$pop_raw   <- as.numeric(df[[COL_VALUE]])

# Filter to ages 25+ and map to project bins
df <- df[!is.na(df$age_lower) & df$age_lower >= 25, ]
df$age <- age_5yr_to_bin[as.character(df$age_lower)]

# Aggregate 5-year groups into project bins
result <- df %>%
  filter(!is.na(age), !is.na(state), !is.na(pop_raw), !is.na(year)) %>%
  group_by(state, year, age) %>%
  summarise(population = sum(pop_raw, na.rm = TRUE), .groups = "drop") %>%
  arrange(state, year, age)

write.csv(result, "data/ibge_population_projections.csv", row.names = FALSE)

cat("\nSaved:", nrow(result), "rows to data/ibge_population_projections.csv\n")
cat("States:", length(unique(result$state)),
    "| Years:", min(result$year), "-", max(result$year),
    "| Age bins:", paste(sort(unique(result$age)), collapse = ", "), "\n")

# Quick sanity check
cat("\nSample (SP, 2030):\n")
print(result %>% filter(state == "SP", year == 2030))
