#' Built-in Manuscript Profiles
#'
#' @return Named list of profile configurations.
mrmoss_profiles <- function() {
  list(
    negative_control = list(
      exposures = c(
        "Ground_coffee_consumption", "Overall_healthy_diet", "Fish_consumption",
        "Fruit_consumption", "Healthy_food_consumption", "Meat_consumption",
        "Psychoactive_drinks_consumption", "Vegetables_consumption", "Beef_consumption",
        "Beer_or_cider_consumption", "Spread_on_bread_consumption", "Bread_consumption",
        "alcohol_consumption_quality", "coffee_consumption_measurement", "Instant_coffee_consumption",
        "Champagne_or_white_wine_consumption", "Cheese_consumption", "Cooked_vegetables_consumption",
        "Decaffeinated_coffee_consumption", "Dried_fruit_consumption", "Drink_temperature",
        "Fresh_fruit_consumption", "Lamb_consumption",
        "Percentage_fat_in_milk_consumption", "Non_oily_fish_consumption", "Oily_fish_consumption",
        "Pork_consumption", "Processed_meat_consumption",
        "Red_wine_consumption", "Salad_consumption", "Added_salt_consumption",
        "Tea_consumption", "Water_consumption_corrected_for_coffee", "Water_consumption"
      ),
      outcomes = c("Emotional_neglect", "Physical_abuse", "Sexual_abuse"),
      iv_thresholds = c(5e-07, 5e-08),
      rd = 1.2,
      output_prefix = "outcome346_D1.2"
    ),
    positive_control = list(
      exposures = c(
        "Apolipoprotein_A_levels",
        "Apolipoprotein_B_levels",
        "HDL_cholesterol_levels",
        "LDL_direct_levels",
        "Lipoprotein_A_levels"
      ),
      outcomes = c(
        "Ischemic_Heart_Disease",
        "Unstable_angina",
        "Myocardial_infarction",
        "Angina_pectoris",
        "Coronary_atherosclerosis"
      ),
      iv_thresholds = c(5e-07, 5e-08),
      rd = 1.2,
      output_prefix = "MVP_outcomes_12345_rd1.2"
    ),
    amd_application = list(
      exposures = c(
        "AgeOfInitiation", "Apolipoprotein_A_levels", "Apolipoprotein_B_levels", "BMI",
        "Body_fat_percentage", "Chronotype", "CigarettesPerDay", "HDL_cholesterol_levels",
        "LDL_direct_levels", "Leisure_screen_time", "Lipoprotein_A_levels", "MVPA",
        "Sedentary_behavior", "Sleep_duration", "SmokingCessation", "Smoking_initiation",
        "Total_cholesterol", "Triglycerides_levels", "WHR", "Walking_pace"
      ),
      outcomes = c("AMD", "AMD_wet", "AMD_dry"),
      iv_thresholds = c(5e-08),
      rd = 1.2,
      output_prefix = "MRMOSS_rd1.2_lifestyle_and_lipids"
    )
  )
}
