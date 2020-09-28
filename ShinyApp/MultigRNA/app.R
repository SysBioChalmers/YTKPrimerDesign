library(shiny)
library(DT)
library(stringr)
library(pracma)


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

reverse = function(s) {
  ss = strsplit(s, "")[[1]]
  revss <- rev(ss)
  return (paste(revss, collapse=""))
}






ui <- fluidPage(
  titlePanel("MultigRNA"), 
  p("Here goes a nice description of the tool"),
  fluidRow(
    column(3,
           textInput("part1", "Sequence part 1", "CAACAAGTTATTACTCAAA", width="100%"),
           textInput("part5", "Sequence part 5", "TTGGTTGATGATAGATTCA", width="100%"),
           textInput("part9", "Sequence part 9", "TTCAAATCTCTGATTCTAT", width="100%")
    ),
    column(3,
           textInput("part2", "Sequence part 2", "AGAAATCTGAAGGTTTGGC", width="100%"),
           textInput("part6", "Sequence part 6", "CTACTGATTCTACAAAAAG", width="100%"),
           textInput("part10", "Sequence part 10", "ATCACAAGAACAACCAATC", width="100%")
    ),
    column(3,
           textInput("part3", "Sequence part 3", "GATAAGGATCCAATTAAAG", width="100%"),
           textInput("part7", "Sequence part 7", "GCAACTGCTGCAATGGCAA", width="100%"),
           textInput("part11", "Sequence part 11", "ACACCATTAGTTAAGGAAC", width="100%"),
           
    ),
    column(3,
           textInput("part4", "Sequence part 4", "CAAGATACATGGTATTACC", width="100%"),
           textInput("part8", "Sequence part 8", "CACCAGATATCGCTAACGA", width="100%"),
           textInput("part12", "Sequence part 12", "TGCAATCAATTAATTTGCC", width="100%")
    )
  ),
  fluidRow(
    column(3,
           textInput("stickyStart", "Start sticky end sequence", "ttta", width="100%")
    ),
    column(3,
           textInput("stickyEnd", "End sticky end sequence", "tggc", width="100%")
    ),
    column(3,
           textInput("sideMargin", "Min. part split size", "3", width="100%")
    ),
    column(3)
  ),
  fluidRow(
    column(6,
           actionButton("genButton", "Generate primers",  width = "100%", height="100%")
    )
  ),
  fluidRow(
    column(12, DT::dataTableOutput("outpPrim", width="100%"))
  )
)

#We place the data in global variables to speed up the recursive function
partStrings = NULL
startPart = NULL
endPart = NULL
cutSeqs = NULL
pieces = NULL
conflPieces = NULL
numParts = NULL
startTag = NULL
endTag = NULL
sideMargin = NULL

#reads data from gui and writes to global vars
setGlobalVars = function(input) {
  partStringsLoc = toupper(
    c(input$part1,
      input$part2,
      input$part3,
      input$part4,
      input$part5,
      input$part6,
      input$part7,
      input$part8,
      input$part9,
      input$part10,
      input$part11,
      input$part12))
  partStringsLoc = partStringsLoc[partStringsLoc != ""]
  if (length(partStringsLoc) <= 1) {
    return (list(FALSE, "You need to add at least 2 sequences"))
  }
  
  
  #so, the first and last should not be split, so they can be removed and stored separately:
  startPart <<- partStringsLoc[1]
  endPart <<- partStringsLoc[length(partStringsLoc)]
  partStrings <<- partStringsLoc[c(-1, -length(partStringsLoc))]
  numParts <<- length(partStrings)
  
  
  startTag <<- toupper(input$stickyStart)
  endTag <<- toupper(input$stickyEnd)
  sideMargin <<- as.numeric(input$sideMargin)
  
  return (list(TRUE, ""))
}


prepareFindSeqCombData = function() {
  #prepare global variables
  #########################
  if (length(partStrings) != 0)
  {
    #first create all possible sequences and their conflicting pieces of length 3
    cutSeqsLoc = vector(mode = "list", length = length(partStrings)) #the 4-letter cut site
    piecesLoc = vector(mode = "list", length = length(partStrings)) #the 2 three letter pieces
    conflPiecesLoc = vector(mode = "list", length = length(partStrings)) #the 4 three letter pieces that needs to be checked against the pieces of other cutSeqs
    for (part in 1:length(partStrings)) {
      partLen = str_length(partStrings[part])
      sm = min(sideMargin, floor((partLen -4)/2)) #make sure at least 1 index is available
      print(partLen)
      print(sm)
      numCuts = partLen-3-sm*2 #don't allow cuts at the edges, need at least 4 nucleotides there
      cs = rep("", numCuts)
      pcs = vector(mode = "list", length = numCuts)
      confl = vector(mode = "list", length = numCuts)
      
      #sort the cuts in a better order, starting from the middle
      center = ceiling((partLen-3)/2)
      indices = center
      for (addOn in 1:ceiling((partLen-3)/2)) { #may loop outside seq, but that is ok as long as it includes all
        num = center + addOn
        if (num <= partLen - 3 - sm) {
          indices = c(indices,num)
        }
        num = center - addOn
        if (num > sm) {
          indices = c(indices,num)
        }
      }
      #print(sm)
      #print(indices)
      #print(numCuts)
      
      for (j in 1:numCuts) {
        i = indices[j]
        cs[[j]] = substr(partStrings[part], i, i+3)
        pcs[[j]] = c(substr(partStrings[part], i, i+2), substr(partStrings[part], i+1, i+3))
        confl[[j]] = c(substr(partStrings[part], i, i+2), substr(partStrings[part], i+1, i+3), 
                       complementary(reverse(substr(partStrings[part], i, i+2))), complementary(reverse(substr(partStrings[part], i+1, i+3))))
      }
      cutSeqsLoc[[part]] = cs
      piecesLoc[[part]] = pcs
      conflPiecesLoc[[part]] = confl
    }
    
    #so, I realize this is a bit odd, but practical in this case... (write to global vars)
    cutSeqs <<- cutSeqsLoc
    pieces <<- piecesLoc
    conflPieces <<- conflPiecesLoc
  } else {
    cutSeqs <<- NULL
    pieces <<- NULL
    conflPieces <<- NULL
  }
}

filterFindSeqCombData = function() {
  if (length(partStrings) != 0)
  {
    #remove cutSeqs that are not compatible with themselves (such as AATT) or that conflicts with the start and end tags
    startEndPieces = c(substr(startTag, 1, 3), substr(startTag, 2, 4), substr(endTag, 1, 3), substr(endTag, 2, 4) )
    
    cutSeqsLoc = cutSeqs
    piecesLoc = pieces
    conflPiecesLoc = conflPieces
    
    for (part in 1:length(partStrings)) {
      toRem = NULL
      numCuts = length(cutSeqsLoc[[part]])
      pcs = piecesLoc[[part]]
      confl = conflPiecesLoc[[part]]
      for (i in 1:numCuts) {
        if (any(pcs[[i]] %in% confl[[i]][3:4]) | any(startEndPieces %in% confl[[i]][3:4])) {
          toRem = c(toRem, i)
          print(paste0("Removing: ", cutSeqsLoc[[part]][[i]]))
        } 
      }
      if (length(toRem) > 0) {
        cutSeqsLoc[[part]] = cutSeqsLoc[[part]][-toRem]
        piecesLoc[[part]] = piecesLoc[[part]][-toRem]
        conflPiecesLoc[[part]] = conflPiecesLoc[[part]][-toRem]
      }
    }
    
    #now write to global vars
    cutSeqs <<- cutSeqsLoc
    pieces <<- piecesLoc
    conflPieces <<- conflPiecesLoc
  }
}

#The recursive function to find the right combination of places to cut the parts
findSeqComb = function(selVector = NULL) {
  currPart = length(selVector) + 1
  #  if (currPart == 6) {
  #    print("test")
  #  }
  if (currPart > numParts) {
    return (list(T, selVector)) #we have succeeded, there are no more to find
  }
  confl = conflPieces[[currPart]]
  numCuts = length(cutSeqs[[currPart]])
  for(i in 1:numCuts) {
    #check if this cut is compatible with the previous
    success = T
    if (!is.null(selVector)) {
      for (j in 1:length(selVector)) {
#        print (selVector)
#        print(j)
        success = success & !any(pieces[[j]][[selVector[j]]] %in% confl[[i]])
      }
    }
    if (success) {
      x = findSeqComb(c(selVector, i))
      if (x[[1]]) {
        return (x) #success
      }
    }
  }
  
  #none found, failure
  return (list(F))
}


extractSeqs = function(algRes) {
  if (length(partStrings) != 0)
  {
    seqs = rep("", length(partStrings))
    for (i in 1:length(partStrings)) {
      seqs[i] = cutSeqs[[i]][[ algRes[[2]] [i] ]]
    }
    
    return(seqs)
  } else {
    return (NULL)
  }
}

extractPositions = function(seqs) {
  if (length(partStrings) != 0)
  {
    positions = rep(NA, length(partStrings))
    for (i in 1:length(partStrings)) {
      res = strfind(partStrings[[i]], seqs[[i]])
      print(partStrings[[i]])
      print(seqs[[i]])
      print(res)
      #take the pos closest to the middle (a sequence could appear several times in a part)
      best = NA
      bestDistFromMiddle = 10000
      middle = ceiling(str_length(partStrings[[i]])/2)
      for (pos in res) {
        dist = abs(pos-middle)
        if (dist < bestDistFromMiddle) {
          bestDistFromMiddle = dist
          best = pos
        }
      }
      positions[i] = best
    }
    
    return(positions)
  } else {
    return (NULL)
  }
}

fwdFill = "gcaf" #just space filler
revFill = "atgc" #just space filler
enz = "cgtctc" # recognition site for the restriction enzyme
standardFwdEnd = "GTTTTAGAGCTAGAAATAGCAAGTTAA"
standardRevEnd = "TTTCTTAGCTGCCTATACGGCAG"

genPrimers = function(positions) {
  
  #so, the first fragment will contain the whole first part (in a separate variable) and half of the first in the list (i+1)
  #the middle fragments will contain the second half of the i part and the first half of the i+1 part
  #the last fragment will contain the first half of the i part and the whole end part (in a separate variable)
  #the first fragment should also contain the start tag sticky end and the last the end tag
  fwdPrimers = NULL
  revPrimers = NULL
  
  for (i in 0:length(partStrings)) {
    if (i == 0) {
      leftPart = startPart
      stt = paste0(tolower(startTag), "a")
    } else {
      leftPart = substr(partStrings[[i]], positions[[i]], str_length(partStrings[[i]]))
      print(leftPart)
      print(positions[[i]])
      print(length(partStrings[[i]]))
      stt = ""
    }
    if (i == length(partStrings)) {
      rightPart = endPart
      et = paste0(tolower(endTag), "a")
    } else {
      rightPart = substr(partStrings[[i+1]], 1, positions[[i+1]] + 3)
      et = ""
    }
    fwd = paste0(fwdFill, enz, "a", stt, "<b>", tolower(leftPart), "</b>", standardFwdEnd)
    rev = paste0(revFill, enz, "a", et, "<b>", tolower(reverse(complementary(toupper(rightPart)))), "</b>", standardRevEnd)
    
    fwdPrimers = c(fwdPrimers, fwd)
    revPrimers = c(revPrimers, rev)
  }
  
  return(list(fwdPrimers, revPrimers))
  
  #first forward
  #fwd1 = paste0(fwdFill, enz, "a", toLower(startTag), substring(partStrings[[1]], 1, positions[[1]]+3), standardFwdEnd)
  #rev1 = paste0(revFill, enz, "a", toLower(startTag), substring(partStrings[[1]], 1, positions[[1]]+3), standardRevEnd)
}


testStuff = function() {
  
  partStrings <<- c("ACGTTCAATTCCCGAAGCGG", "AAGAATTGTTCCC")
  numParts <<- length(partStrings)
  startTag <<- "CCAA"
  endTag <<- "GTTG"
  sideMargin <<- 3
  cutSeqs <<- NULL
  prepareFindSeqCombData()
  
  expRes = c("TTCC", "TCCC", "ATTC", "CCCG", "AATT", "CCGA", "CAAT", "CGAA", "TCAA", "GAAG", "TTCA")
  expRes2 = c("ATTG", "TTGT", "AATT", "TGTT")
  #print(cutSeqs[[1]])
  #print(cutSeqs[[2]])
  #check cutSeqs
  print(all(cutSeqs[[1]] == expRes))
  print(all(cutSeqs[[2]] == expRes2))
  #check pieces (only check 1 item)
  all(pieces[[1]][[6]] == c("CCG", "CGA"))
  #check conflPieces (only check 1 item)
  all(conflPieces[[1]][[6]] == c("CCG", "CGA", "CGG", "TCG"))
  
  #now test filtering
  filterFindSeqCombData()
  print(cutSeqs[[1]])
  expRes = c("TTCC", "TCCC", "ATTC", "CCCG", "CCGA", "CGAA", "GAAG", "TTCA")
  expRes2 = c("TGTT")
  print(all(cutSeqs[[1]] == expRes))
  print(all(cutSeqs[[2]] == expRes2))
  
  
  #construct a minimal case with only one solution:
  #partStrings <<- c("ACGATCCCCGG", "AAGCCCTGCCC", "AAGCCCTACCACCC", "AAGATCCCCC")
  partStrings <<- c("ACGATCCCCGG", "AAGTGCCCCCC", "AAGCCCTACCACCC", "AAGCTACCCC")
  numParts <<- length(partStrings)
  prepareFindSeqCombData()
  res = findSeqComb()
  seqs = extractSeqs(res)
  #we just check that the results it spits out is a solution without conflicts, 
  #which the following is:
  expRes = c("ATCC", "TGCC", "ACCA", "CTAC")
  all(seqs == expRes)
  
  #now find where this is in the parts:
  positions = extractPositions(seqs)
  #test by generating the sequences from the positions and compare to expected result:
  newSeqs = rep(NA,length(positions))
  for (i in 1:length(positions)) {
    newSeqs[[i]] = substr(partStrings[[i]], positions[i], positions[i]+3)
  }
  all(expRes == newSeqs)
  
  
}


onGenButton = function(input, output, session) {
  res = setGlobalVars(input)
  if (res[[1]]) {
    prepareFindSeqCombData()
    filterFindSeqCombData()
    
    # All global variables set, run the recursive function to find the place to cut each part
    res = findSeqComb()
    if (res[[1]]) {
    
      #extract the resulting sticky end sequences
      seqs = extractSeqs(res)
  
      #now find where this is in the parts:
      positions = extractPositions(seqs)
      print(positions)
      
      #generate primers:
      primers = genPrimers(positions)
      dfPrim = data.frame(Forward = primers[[1]], Reverse = primers[[2]], stringsAsFactors = F)
      colnames(dfPrim) = c("Forward primer", "Reverse primer")
      
      output$outpPrim = DT::renderDataTable(DT::datatable({dfPrim}, escape = F, width="100%", rownames = FALSE, class="compact", selection = 'none', options = list(dom = 't', pageLength = 20)))
    } else {
      showModal(modalDialog(
        title = "Generation failure.",
        "Couldn't find unique internal sticky ends.",
        easyClose = TRUE,
        footer = NULL
      ))
      output$outpPrim = DT::renderDataTable(DT::datatable({NULL}, escape = F, width="100%", rownames = FALSE, class="compact", selection = 'none', options = list(dom = 't', pageLength = 20)))
      
    }
  } else {
    showModal(modalDialog(
      title = "Generation failure.",
      res[[2]],
      easyClose = TRUE,
      footer = NULL
    ))
    output$outpPrim = DT::renderDataTable(DT::datatable({NULL}, escape = F, width="100%", rownames = FALSE, class="compact", selection = 'none', options = list(dom = 't', pageLength = 20)))
  }
}



server <- function(input, output, session) {
  observeEvent(input$genButton, {
    onGenButton(input, output, session)
  })
  
}



#don't add any code after the shinyApp line below, it has to be the last line!
shinyApp(ui = ui, server = server)
