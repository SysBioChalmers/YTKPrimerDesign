#install.packages('rsconnect')
#the line below needs to be run once, and will save some token files to disk
#rsconnect::setAccountInfo(name='johan-gson', token='8F82E4588E476DE913B8A5A098394D9D', secret='<fill in the secret, can be collected from the shinyapps.io>')

setwd("C:/Work/MatlabCode/projects/YTKPrimerDesign/YTKPrimerDesign/ShinyApp/The other tool")
rsconnect::setAccountInfo(name='johan-gson', token='8F82E4588E476DE913B8A5A098394D9D', secret='<fill in>')


library(rsconnect)
rsconnect::deployApp('C:/Work/MatlabCode/projects/YTKPrimerDesign/YTKPrimerDesign/ShinyApp/YTK Primer Design')
rsconnect::deployApp('C:/Work/MatlabCode/projects/YTKPrimerDesign/YTKPrimerDesign/ShinyApp/MultigRNA', appTitle = 'MultigRNA')
