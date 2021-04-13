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
  if (s == "") {
    return ("")
  }
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

startValueNames = c("Type1", "Type 2", "Type 3", "Type 3A", "Type 3B", "Type 4", "Type 4A", "Type 4B", "Type 5", "Type 6", "Type 7", "Type 8", "Type 8A", "Type 8B")
startValueVals = c("CCCT", "AACG", "TATG", "TATG", "TTCT", "ATCC", "ATCC", "TGGC", "GCTG", "TACA", "GAGT", "CCGA", "CCGA", "CAAT")
startValueExtras = c("GAATTCGCATCTAGA", "", "", "", "", "CTCGAG", "", "", "", "", "", "GCGGCCGC", "GCGGCCGC", "")
startValueStopCodons = c("", "", "", "", "", "TAA", "", "", "", "", "", "", "", "")


endValueNames = c("Type 1", "Type 2", "Type 3", "Type 3A", "Type 3B", "Type 4", "Type 4A", "Type 4B", "Type 5", "Type 6", "Type 7", "Type 8", "Type 8A", "Type 8B")
endValueVals = c("AACG", "TATG", "ATCC", "TTCT", "ATCC", "GCTG", "TGGC", "GCTG", "TACA", "GAGT", "CCGA", "CCCT", "CAAT", "CCCT")
endValueExtras = c("", "AGATC", "GG", "GG", "GG", "", "CTCGAG", "", "ACTAGTGCACTGCAG",  "", "", "GCGGCCGC", "GCGGCCGC", "")
endValueStopCodons = c("", "", "", "", "", "", "TAA", "", "", "", "", "", "", "")

enzymeRecSitesIds = c("BsaI", "BsmBI", "NotI", "EcoRI", "XbaI", "SpeI", "PstI", "BglII", "BamHI", "XhoI")
enzymeRecSitesSeqs = c("GGTCTC", "CGTCTC", "GCGGCCGC", "GAATTC", "TCTAGA", "ACTAGT", "CTGCAG", "AGATCT", "GGATCC", "CTCGAG")
enzymeRecExtraTexts = c("", "", " (NotI is used for linearisation of integration cassettes)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)",
                        " (Restriction site should be removed if part will also be used as BioBrick or BglBrick)")
                        

#startChoices = startValueVals
#names(startChoices) = startValueNames
#endChoices = endValueVals
#names(endChoices) = endValueNames




ui <- fluidPage(
  titlePanel("Yeast Toolkit Primer Design Tool"), 
  tabsetPanel(
    id = 'mainTab',
    tabPanel("Application", 
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
      )
    ),
    tabPanel("Instructions", 
            h2("Summary"),
            p(HTML("This tool was published in 'ourPaperReference'. It can be used to design new genetic parts for modular cloning in <i>Saccharomyces cerevisiae</i>. The underlying methodology was developed by Lee et al. (2015). The design tool will provide a list of primer sequences with type-specific overhangs. Additionally, the sequence of the Golden Gate Assembly (Engler et al. 2008) of the PCR product and the entry vector pYTK001 can be downloaded as a FASTA file. For more information see Lee et al. (2015) and 'ourPaperReference'.<br>")),
            h2("Instructions"),
            p(HTML("  1. Paste the DNA sequence you want to include in the toolkit<br>")),
            p(HTML("  2. Choose the desired type-specific overhangs from the dropdown menus<br>")),
            p(HTML("  3. Press 'Gen. primers' to generate a list of forward and reverse primers. The user can pick a primer pair depending on the length, GC-content and melting temperature of the template-binding sequence<br>")),
            p(HTML("  4. A warning will be displayed if BsaI, BsmBI, EcoRI, XbaI, SpeI, PstI, BglII, BamHI, XhoI or NotI recognition sites are found<br>")),
            p(HTML("  5. Download the plasmid sequence (Golden Gate Assembly of PCR product and pYTK001 using BsmBI) by pressing 'Export FASTA'<br>")),
            h2("References"),
            p(HTML("Our reference<br>")),
            p(HTML('Lee, M. E., et al. (2015). "A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly." ACS Synthetic Biology 4(9): 975-986.<br>')),
            p(HTML('Engler, C., et al. (2008). "A One Pot, One Step, Precision Cloning Method with High Throughput Capability." PLoS One 3(11): e3647.<br>')),
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
      revPrimTot = paste0(tolower(complementary(reverse(univEndTag))), "<b>", tolower(reverse(complementary(endTag))), "</b>", tolower(reverse(complementary(endValueExtra))), tolower(reverse(complementary(endValueStopCodon))), complementary(reversePrimer)) 
      dfRev = rbind(dfRev, c(revPrimTot,
                             round(gcContent(reversePrimer), 2),
                             le,
                             round(meltTemp(reversePrimer), 1)), stringsAsFactors=F)
    }
  }
#  print(forwardPrimer)
#  print(reversePrimer)
  
  colnames(dfFwd) = c("Primer sequence", "GC", "Len", "Tm")
  colnames(dfRev) = c("Primer sequence", "GC", "Len", "Tm")
  
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
                    "Start codon is already included in fwd primer. Do not include stop codon in pasted sequence if part will be used as fusion or tagged protein. (see Lee et al. 2015)\n")
  }
  if (input$endCombo == "Type 3B") {
    strWarn = paste0(strWarn,
                    "Do not include stop codon in pasted sequence if part will be used as fusion or tagged protein. (see Lee et al. 2015)\n")
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
    filename = function(){ paste0(input$filename, ".FASTA")
    },
    #filename = paste0("test", ".FASTA"),
    content = function(file) {
      stringsToWrite = paste0("> ", input$startCombo, "-", input$filename, "-", input$endCombo)
      
      #generate sequence
      fwdOverhangStd = "TCGGTCTCA"
      fwdOverhangTypeSpec = lookup(input$startCombo, startValueNames, startValueVals)
      fwdOverhangTypeSpecExtra = lookup(input$startCombo, startValueNames, startValueExtras)
      fwdOverhangTypeStopCodon = lookup(input$startCombo, startValueNames, startValueStopCodons)
      seq = cleanupSeq(input$inputSeq)
      revOverhangStd = "TGAGACC"
      revOverhangTypeSpec = lookup(input$endCombo, endValueNames, endValueVals)
      revOverhangTypeSpecExtra = lookup(input$endCombo, endValueNames, endValueExtras)
      revOverhangTypeStopCodon = lookup(input$endCombo, endValueNames, endValueStopCodons)
      
      backbone = "AGACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCAGATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCCAAGCGAGCTCGATATCAAATTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAAACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGAAATTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCCGTAATATCCAGCTGAACGGTCTGGTTATAGGTACATTGAGCAACTGACTGAAATGCCTCAAAATGTTCTTTACGATGCCATTGGGATATATCAACGGTGGTATATCCAGTGATTTTTTTCTCCATTTTAGCTTCCTTAGCTCCTGAAAATCTCGATAACTCAAAAAATACGCCCGGTAGTGATCTTATTTCATTATGGTGAAAGTTGGAACCTCTTACGTGCCCGATCAATCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTAG"
      
      backbone8_8A = "AGTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATACTGAGCACGGATCTGAAAGAGGAGAAAGGATCTATGGCGAGTAGCGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACTCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAAGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACTAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAAGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAAGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAATAAGGATCTCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATAAG"
      #if Type 8 or 8A is selected for either prefix or suffix, another backbone should be selected:
      if (input$startCombo == "Type 8" | input$startCombo == "Type 8A" | input$endCombo == "Type 8" | input$endCombo == "Type 8A") {
        backbone = backbone8_8A
      }
      
      totSeqToExp = paste0(fwdOverhangStd, fwdOverhangTypeSpec, fwdOverhangTypeStopCodon, fwdOverhangTypeSpecExtra, seq, revOverhangTypeStopCodon, revOverhangTypeSpecExtra, revOverhangTypeSpec, revOverhangStd, backbone)

      
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
