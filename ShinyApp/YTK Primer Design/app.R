library(shiny)
library(DT)
library(stringr)
#install.packages("shinyBS")
library(shinyBS)
#install.packages("pracma")
library(pracma)
library(qdapTools)

#runExample("01_hello")

#there's a bug in textAreaInput, fixed by the code below (which is a copy of the same input)
#See https://stackoverflow.com/questions/40808889/a-dynamically-resizing-shiny-textareainput-box
#based on Shiny textAreaInput
textAreaInput2 <- function (inputId, label, value = "", width = NULL, height = NULL, 
                            cols = NULL, rows = NULL, placeholder = NULL, resize = NULL) 
{
  value <- restoreInput(id = inputId, default = value)
  if (!is.null(resize)) {
    resize <- match.arg(resize, c("both", "none", "vertical", 
                                  "horizontal"))
  }
  style <- paste("max-width: 100%;", if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"), if (!is.null(height)) 
      paste0("height: ", validateCssUnit(height), ";"), if (!is.null(resize)) 
        paste0("resize: ", resize, ";"))
  if (length(style) == 0) 
    style <- NULL
  div(class = "form-group", 
      tags$label(label, `for` = inputId), tags$textarea(id = inputId, 
                                                        class = "form-control", placeholder = placeholder, style = style, 
                                                        rows = rows, cols = cols, value))
}

#super inefficent, but that doesn't really matter
cleanupSeq = function(s) {
  s = toupper(s)
  if (str_length(s) == 0) {
    return(s)
  }
  ss = strsplit(s,"")[[1]]
  ret = ""
  for (i in 1:length(ss)) { 
    if ((ss[i] == 'A') | (ss[i] == 'T') | (ss[i] == 'C') | (ss[i] == 'G')) {
      ret = paste0(ret,ss[i])
    }
  }
  return (ret)
}

compChar = function(s) {
  if (s == "A") return ("T")
  if (s == "T") return ("A")
  if (s == "C") return ("G")
  if (s == "G") return ("C")
  return("X")
}


complementary = function(s) {
  ss = strsplit(s,"")[[1]]
  ret = ""
  for (i in 1:length(ss)) {
    ret = paste0(ret, compChar(ss[i]))
  }
  return (ret)
}

gcContent = function(s) {
  return(str_count(s, "[GC]")/str_length(s))
}

#Wallace formula, 
#Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
#Wallace formula: Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
#Wallace RB et al. (1979) Nucleic Acids Res 6:3543-3557, PMID 158748
meltTemp = function(s) {
  wA = str_count(s, "A")
  xT = str_count(s, "T")
  yG = str_count(s, "G")
  zC = str_count(s, "C")

  return (64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC))
}

checkPrimer = function(s) {
  mt = meltTemp(s)
  return ((mt >= 48) & (mt <= 60))
}

reverse = function(s) {
  ss = strsplit(s, "")[[1]]
  revss <- rev(ss)
  return (paste(revss, collapse=""))
}


univStartTag = "GCATCGTCTCATCGGTCTCA"
univEndTag = "TGAGACCTGAGACGGCAT"

startValueNames = c("Type 2", "Type 3", "Type 3A", "Type 3B", "Type 4", "Type 4A", "Type 4B", "Type 6", "Type 7", "Type 8", "Type 8A", "Type 8B")
startValueVals = c("AACG", "TATG", "TATG", "TTCT", "ATCC", "ATCC", "TGGC", "TACA", "GAGT", "CCGA", "CCGA", "CAAT")
startValueExtras = c("", "", "", "", "CTCGAG", "", "", "", "", "GCGGCCGC", "GCGGCCGC", "")
startValueStopCodons = c("", "", "", "", "TAA", "", "", "", "", "", "", "")


endValueNames = c("Type 2", "Type 3", "Type 3A", "Type 3B", "Type 4", "Type 4A", "Type 4B", "Type 6", "Type 7", "Type 8", "Type 8A", "Type 8B")
endValueVals = c("TATG", "ATCC", "TTCT", "ATCC", "GCTG", "TGGC", "GCTG", "GAGT", "CCGA", "CCCT", "CAAT", "CCCT")
endValueExtras = c("AGATC", "GG", "GG", "GG", "", "CTCGAG", "", "", "", "GCGGCCGC", "GCGGCCGC", "")
endValueStopCodons = c("", "", "", "", "", "TAA", "", "", "", "", "", "")

enzymeRecSitesIds = c("BsaI", "BsmBI", "NotI", "BglII", "BamHI", "XhoI")
enzymeRecSitesSeqs = c("GGTCTC", "CGTCTC", "GCGGCCGC", "AGATCT", "GGATCC", "CTCGAG")
enzymeRecExtraTexts = c("", "", " (NotI is used for linearisation of integration cassettes)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)")
                        

#startChoices = startValueVals
#names(startChoices) = startValueNames
#endChoices = endValueVals
#names(endChoices) = endValueNames




ui <- fluidPage(
  titlePanel("Yeast Toolkit Primer Design Tool"), 
  p("Here goes a nice description of the tool"),
  fluidRow(
    column(12,
           textAreaInput2("inputSeq", "Input sequence", width="100%", resize="none", cols=140, rows=8)
    )
  ),
  fluidRow(
    column(6,
           mainPanel(width="100%",
                     fluidRow(
                       column(6, selectInput("startCombo", "Part type prefix", startValueNames,  width = "100%", selectize = F)),
                       column(6, selectInput("endCombo", "Part type suffix", endValueNames,  width = "100%", selectize = F))
                     ),
                     fluidRow(
                       column(4, checkboxInput("filterCheck", "Filter primers", T)),
                       column(4, actionButton("genButton", "Gen. primers",  width = "100%")),
                       column(4,actionButton("expButton", "Export FASTA",  width = "100%"))
                     )
                     
           )
           
    ),
    column(6,
           uiOutput("warningsTitle"),
           verbatimTextOutput("warnings"),
           tags$head(tags$style("#warnings{color: red;overflow-y:scroll;max-height: 100px;}"))

    )
  ),
  fluidRow(
    column(12,
           conditionalPanel(
             'output.primersGenerated === "t"',
             tabsetPanel(
                id = 'outputTab',
                tabPanel("Forward primer", DT::dataTableOutput("outpFwd", width="100%")),
                tabPanel("Reverse primer", DT::dataTableOutput("outpRev", width="100%"))
             )
           )
    )
  ),
  bsModal("expModalDlg", "Export FASTA", "expButton", size = "large",textInput('filename', 'Filename', "Sequence"), downloadButton('downloadSeq', 'Download'))
)

primGen = "f"

onGenButton = function(input, output, session) {
  primGen = "t"
  output$primersGenerated = renderText(primGen)
  
  inp = paste0(input$inputSeq)
  sq = cleanupSeq(inp)
  startTag = lookup(input$startCombo, startValueNames, startValueVals)
  startValueExtra = lookup(input$startCombo, startValueNames, startValueExtras)
  startValueStopCodon = lookup(input$startCombo, startValueNames, startValueStopCodons)
  endTag = lookup(input$endCombo, endValueNames, endValueVals)
  endValueExtra = lookup(input$endCombo, endValueNames, endValueExtras)
  endValueStopCodon = lookup(input$endCombo, endValueNames, endValueStopCodons)
  

  dfFwd <- data.frame(Sequence=character(),
                      GCCont=double(), 
                      Length=integer(), 
                      MeltTm=double(),
                      stringsAsFactors=FALSE)   
  dfRev <- data.frame(Sequence=character(),
                      GCCont=double(), 
                      Length=integer(), 
                      MeltTm=double(),
                      stringsAsFactors=FALSE)   
  
  
  for (le in 17:28) {
    forwardPrimer = substr(sq, 1, le)
    if ((!input$filterCheck) | checkPrimer(forwardPrimer)) {
#      fwdPrimTot = paste0(tolower(univStartTag), "<b>", tolower(startTag), '</b><span style="color:red">', tolower(startValueStopCodon), "</span><u>", tolower(startValueExtra), "</u>", forwardPrimer)
      fwdPrimTot = paste0(tolower(univStartTag), "<b>", tolower(startTag), '</b>', tolower(startValueStopCodon), tolower(startValueExtra), forwardPrimer)
      dfFwd = rbind(dfFwd, c(fwdPrimTot,
                             round(gcContent(forwardPrimer), 2),
                             le,
                             round(meltTemp(forwardPrimer), 1)), stringsAsFactors=F)
    }
    
    reversePrimer = reverse(substr(sq, str_length(sq) - le + 1, str_length(sq))) #reverses automatically
    if ((!input$filterCheck) | checkPrimer(reversePrimer)) {
      #revPrimTot = paste0(tolower(complementary(reverse(univEndTag))), '<span style="color:red">', tolower(endValueStopCodon),"</span><u>", tolower(endValueExtra),"</u><b>", tolower(endTag), "</b>", complementary(reversePrimer)) 
      revPrimTot = paste0(tolower(complementary(reverse(univEndTag))), tolower(endValueStopCodon), tolower(endValueExtra),"<b>", tolower(endTag), "</b>", complementary(reversePrimer)) 
      dfRev = rbind(dfRev, c(revPrimTot,
                             round(gcContent(reversePrimer), 2),
                             le,
                             round(meltTemp(reversePrimer), 1)), stringsAsFactors=F)
    }
  }
#  print(forwardPrimer)
#  print(reversePrimer)
  
  colnames(dfFwd) = c("Primer sequence", "GC", "Len", "Melt Tm")
  colnames(dfRev) = c("Primer sequence", "GC", "Len", "Melt Tm")
  
  output$outpFwd = DT::renderDataTable(DT::datatable({dfFwd}, escape = F, width="100%", rownames = FALSE, class="compact", selection = 'none', options = list(dom = 't', pageLength = 20)))
  output$outpRev = DT::renderDataTable(DT::datatable({dfRev}, escape = F, width="100%", rownames = FALSE, class="compact", selection = 'none', options = list(dom = 't', pageLength = 20)))

  #Warnings:
  strWarn = ""
  for (i in 1:length(enzymeRecSitesIds)) {
    ind = strfind(sq, enzymeRecSitesSeqs[i])
    for (a in ind) {
      strWarn = paste0(strWarn, "Restr. enzyme recogn. site found in seq: ", enzymeRecSitesIds[i], ": ", enzymeRecSitesSeqs[i], " at pos ", a, enzymeRecExtraTexts[i], "\n")
    }
    #if the enzymes are not the same after reverse and complement (which some are), look for the reverse complement as well
    #This was tested manually
    revCompl = reverse(complementary(enzymeRecSitesSeqs[i]))
    if (enzymeRecSitesSeqs[i] != revCompl) {
      ind = strfind(sq, revCompl)
      for (a in ind) {
        strWarn = paste0(strWarn, "Restr. enzyme recogn. site found in the reverse compl. seq: ", enzymeRecSitesIds[i], ": ", revCompl, " at pos ", a, enzymeRecExtraTexts[i], "\n")
      }
    }
  }
  #issue warnings for certain parts, since the start/end codon is included in the primer
  if (input$startCombo == "Type 3" | input$startCombo == "Type 3A") {
    strWarn = paste0(strWarn,
                    "Start codon is already included in fwd primer. Do not include stop codon in pasted sequence if part will be used as fusion or tagged protein.\n")
  }
  if (input$endCombo == "Type 3B") {
    strWarn = paste0(strWarn,
                    "Stop codon is already included in rev primer. Do not include stop codon in pasted sequence if part will be used as fusion or tagged protein.\n")
  }
#Test scrollbar
#  for (i in 1:10) {
#    strWarn = paste0(strWarn, "Warning: Test\n")
#  }
  
  #updateTextAreaInput(session, "warnings", value = strWarn)
  output$warningsTitle = renderText({HTML("<b>Warnings</b>")})
  output$warnings = renderText(strWarn)
  
}



server <- function(input, output, session) {
  output$primersGenerated = renderText(primGen)
  outputOptions(output, "primersGenerated", suspendWhenHidden = FALSE)
  observeEvent(input$genButton, {
    onGenButton(input, output, session)
  })

  output$downloadSeq <- downloadHandler(
    filename = paste0(input$filename, ".FASTA"),
    #filename = paste0("test", ".FASTA"),
    content = function(file) {
      stringsToWrite = paste0(">Yeast toolkit export. Name: ", input$filename, " Starts at ", input$startCombo, ", ends at ", input$endCombo)
      
      #generate sequence
      fwdOverhangStd = "TCGGTCTCA"
      fwdOverhangTypeSpec = input$startCombo
      seq = cleanupSeq(input$inputSeq)
      revOverhangStd = "TGAGACC"
      revOverhangTypeSpec = input$endCombo
      backbone = "AGACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCAGATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCCAAGCGAGCTCGATATCAAATTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAAACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGAAATTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCCGTAATATCCAGCTGAACGGTCTGGTTATAGGTACATTGAGCAACTGACTGAAATGCCTCAAAATGTTCTTTACGATGCCATTGGGATATATCAACGGTGGTATATCCAGTGATTTTTTTCTCCATTTTAGCTTCCTTAGCTCCTGAAAATCTCGATAACTCAAAAAATACGCCCGGTAGTGATCTTATTTCATTATGGTGAAAGTTGGAACCTCTTACGTGCCCGATCAATCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTAG"
      
      totSeqToExp = paste0(fwdOverhangStd, fwdOverhangTypeSpec, seq, revOverhangTypeSpec, revOverhangStd, backbone)

      
      ll = str_length(totSeqToExp)
      for (i in seq(1, ll, by=80)) {
        stringsToWrite = c(stringsToWrite,paste0(substr(totSeqToExp, i, i+79)))
      }
      writeLines(stringsToWrite, con=file)
      toggleModal(session, modalId = "expModalDlg") #close the dialog
    }) 
}



#don't add any code after the shinyApp line below, it has to be the last line!
shinyApp(ui = ui, server = server)
