
# Functions that create arrays for cases, deaths and recovered

# Cases
ArrayCases <- function(firstday, df, numbreg) {
 # firstday = first day to include in analysis
 # df = dataframe used (in wide format)
 # numbreg = number of geographical regions
  MondayCases <- df[, seq(firstday, ncol(df), 7)]
  TuesdayCases <- df[, seq((firstday+1), ncol(df), 7)]
  WednesdayCases <- df[, seq((firstday+2), ncol(df), 7)]
  ThursdayCases <- df[, seq((firstday+3), ncol(df), 7)]
  FridayCases <- df[, seq((firstday+4), ncol(df), 7)]
  SaturdayCases <- df[, seq((firstday+5), ncol(df), 7)]
  SundayCases <- df[, seq((firstday+6), ncol(df), 7)]
  
  cases_array <- array(c(unlist(MondayCases), unlist(TuesdayCases), 
                         unlist(WednesdayCases), unlist(ThursdayCases), 
                         unlist(FridayCases), unlist(SaturdayCases),
                         unlist(SundayCases)), c(numbreg, ncol(MondayCases), 7))
  return(cases_array)
}

# Deaths
ArrayDeaths <- function(firstday, df, numbreg) {
  # firstday = first day to include in analysis
  # df = dataframe used (in wide format)
  # numbreg = number of geographical regions
  MondayDeaths <- df[, seq(firstday, ncol(df), 7)]
  TuesdayDeaths <- df[, seq((firstday+1), ncol(df), 7)]
  WednesdayDeaths <- df[, seq((firstday+2), ncol(df), 7)]
  ThursdayDeaths <- df[, seq((firstday+3), ncol(df), 7)]
  FridayDeaths <- df[, seq((firstday+4), ncol(df), 7)]
  SaturdayDeaths <- df[, seq((firstday+5), ncol(df), 7)]
  SundayDeaths <- df[, seq((firstday+6), ncol(df), 7)]

  deaths_array <- array(c(unlist(MondayDeaths), unlist(TuesdayDeaths), 
                         unlist(WednesdayDeaths), unlist(ThursdayDeaths), 
                         unlist(FridayDeaths), unlist(SaturdayDeaths),
                         unlist(SundayDeaths)), c(numbreg, ncol(MondayDeaths), 7))
  return(deaths_array)
}

# Recovered
ArrayRecovered <- function(firstday, df, numbreg) {
  # firstday = first day to include in analysis
  # df = dataframe used (in wide format)
  # numbreg = number of geographical regions
  MondayRec <- df[, seq(firstday, ncol(df), 7)]
  TuesdayRec <- df[, seq((firstday+1), ncol(df), 7)]
  WednesdayRec <- df[, seq((firstday+2), ncol(df), 7)]
  ThursdayRec <- df[, seq((firstday+3), ncol(df), 7)]
  FridayRec <- df[, seq((firstday+4), ncol(df), 7)]
  SaturdayRec <- df[, seq((firstday+5), ncol(df), 7)]
  SundayRec <- df[, seq((firstday+6), ncol(df), 7)]
  
  rec_array <- array(c(unlist(MondayRec), unlist(TuesdayRec), 
                          unlist(WednesdayRec), unlist(ThursdayRec), 
                          unlist(FridayRec), unlist(SaturdayRec),
                          unlist(SundayRec)), c(numbreg, ncol(MondayRec), 7))
  return(rec_array)
}

