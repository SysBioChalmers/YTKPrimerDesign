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


enzymeRecSitesIds = c("BsmBI", "BsaI")
enzymeRecSitesSeqs = c("CGTCTC", "GGTCTC")




ui <- fluidPage(
  titlePanel("MultigRNA"), 
  tabsetPanel(
    id = 'mainTab',
    tabPanel("Application",
      p(HTML("<br>")),
      fluidRow(
        column(6,
               textInput("FwdBP", "Forward primer binding part", "GTTTTAGAGCTAGAAATAGCAAGTTAA", width="100%"),
        ),
        column(6,
               textInput("RevBP", "Reverse primer binding part", "TTTCTTAGCTGCCTATACGGCAG", width="100%"),
        )
      ),      
      fluidRow(
        column(3,
               textInput("part1", "gRNA 1 sequence", width="100%"),
               textInput("part5", "gRNA 5 sequence", width="100%"),
               textInput("part9", "gRNA 9 sequence", width="100%")
        ),
        column(3,
               textInput("part2", "gRNA 2 sequence", width="100%"),
               textInput("part6", "gRNA 6 sequence", width="100%"),
               textInput("part10", "gRNA 10 sequence", width="100%")
        ),
        column(3,
               textInput("part3", "gRNA 3 sequence", width="100%"),
               textInput("part7", "gRNA 7 sequence", width="100%"),
               textInput("part11", "gRNA 11 sequence", width="100%"),
               
        ),
        column(3,
               textInput("part4", "gRNA 4 sequence", width="100%"),
               textInput("part8", "gRNA 8 sequence", width="100%"),
               textInput("part12", "gRNA 12 sequence", width="100%")
        )
    #    fluidRow(
    #      column(3,
    #             textInput("part1", "Sequence part 1", "CAACAAGTTATTACTCAAA", width="100%"),
    #             textInput("part5", "Sequence part 5", "TTGGTTGATGATAGATTCA", width="100%"),
    #             textInput("part9", "Sequence part 9", "TTCAAATCTCTGATTCTAT", width="100%")
    #      ),
    #      column(3,
    #             textInput("part2", "Sequence part 2", "AGAAATCTGAAGGTTTGGC", width="100%"),
    #             textInput("part6", "Sequence part 6", "CTACTGATTCTACAAAAAG", width="100%"),
    #             textInput("part10", "Sequence part 10", "ATCACAAGAACAACCAATC", width="100%")
    #      ),
    #      column(3,
    #             textInput("part3", "Sequence part 3", "GATAAGGATCCAATTAAAG", width="100%"),
    #             textInput("part7", "Sequence part 7", "GCAACTGCTGCAATGGCAA", width="100%"),
    #             textInput("part11", "Sequence part 11", "ACACCATTAGTTAAGGAAC", width="100%"),
    #             
    #      ),
    #      column(3,
    #             textInput("part4", "Sequence part 4", "CAAGATACATGGTATTACC", width="100%"),
    #             textInput("part8", "Sequence part 8", "CACCAGATATCGCTAACGA", width="100%"),
    #             textInput("part12", "Sequence part 12", "TGCAATCAATTAATTTGCC", width="100%")
    #      )
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
        ),
        column(6,
               uiOutput("warningsTitle"),
               verbatimTextOutput("warnings"),
               tags$head(tags$style("#warnings{color: red;overflow-y:scroll;max-height: 100px;}"))
               
        )
      ),
      fluidRow(
        column(12, DT::dataTableOutput("outpPrim", width="100%"))
      )
    ),
    tabPanel("Instructions", 
             h2("Summary"),
             p(HTML("This tool is published in 'our reference'. It is based on the multiplexed gRNA expression with 
                    the use of Csy4 endonuclease (Ferreira et al. 2018). The cloning strategy is based on the vector 
                    pYTK050 from Lee et al. (2015) which is a vector for gRNA cloning via Golden Gate assembly 
                    (Engler et al. 2008). The gRNAs are cloned in the place of a GFP gene or a yeast counterselection 
                    marker. The cloning position has two BsmBI recognition sites and the 3' and 5' sticky ends can vary 
                    according to the destination vector. The present tool designs primers for multiplexed gRNA cloning 
                    as described in 'our reference' and avoids the presence of more than 3 nucleotides in a row in 
                    common in two sticky ends. The template for all PCR reactions is the plasmid pMCL008_28bp which 
                    has a gRNA scaffold followed by a 28bp Csy4 recognition loop. The resulting PCR products can be 
                    used for Golden Gate cloning. This tool can develop primers for cloning up to 12 gRNAs in a single 
                    transcript.<br>")),
             h2("Instructions"),
             p(HTML("  1. Write the starting and ending sticky end sequences that will be created with BsmBI digestion of the destination vector.<br>")),
             p(HTML("  2. Paste the 20bp gRNA sequence in the 'gRNA' boxes. If less than 12 gRNAs will be cloned, leave the rest of the boxes blank.<br>")),
             p(HTML("  3. Click on 'Generate primers'. A pair of primers for each fragment will appear. If any of the gRNAs to be cloned has a BsmBI position a warning will appear.<br>")),
             h2("References"),
             p(HTML("Our reference<br>")),
             p(HTML('Ferreira, Raphael, et al. "Multiplexed CRISPR/Cas9 genome editing and gene regulation using Csy4 in Saccharomyces cerevisiae." ACS synthetic biology 7.1 (2018): 10-15.<br>')),
             p(HTML('Lee, M. E., et al. (2015). "A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly." ACS Synthetic Biology 4(9): 975-986.<br>')),
             p(HTML('Engler, C., et al. (2008). "A One Pot, One Step, Precision Cloning Method with High Throughput Capability." PLoS One 3(11): e3647.<br>')),
    )
  )
)

#We place the data in global variables to speed up the recursive function
partStrings = NULL
startPart = NULL
endPart = NULL
cutSeqs = NULL
pieces = NULL
conflPieces = NULL
pieces2 = NULL
conflPieces2 = NULL
pieces3 = NULL
conflPieces3 = NULL
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
    #conflPiecesLoc = vector(mode = "list", length = length(partStrings)) #the 4 three letter pieces that needs to be checked against the pieces of other cutSeqs
    #So, we should check the following: (1,2,3), (2,3,4), (1,2,4), (1,3,4), and their reverse complements
    #we therefore place (1,2,3), (2,3,4) and their reverse complements in pieces and and conflPieces, and the others the ones named 2 and 3
    piecesLoc = vector(mode = "list", length = length(partStrings)) #the 2 three letter pieces
    conflPiecesLoc = vector(mode = "list", length = length(partStrings)) 
    piecesLoc2 = vector(mode = "list", length = length(partStrings))
    conflPiecesLoc2 = vector(mode = "list", length = length(partStrings)) 
    piecesLoc3 = vector(mode = "list", length = length(partStrings))
    conflPiecesLoc3 = vector(mode = "list", length = length(partStrings)) 
    for (part in 1:length(partStrings)) {
      partLen = str_length(partStrings[part])
      sm = min(sideMargin, floor((partLen -4)/2)) #make sure at least 1 index is available
      print(partLen)
      print(sm)
      numCuts = partLen-3-sm*2 #don't allow cuts at the edges, need at least 4 nucleotides there
      cs = rep("", numCuts)
      pcs = vector(mode = "list", length = numCuts)
      confl = vector(mode = "list", length = numCuts)
      pcs2 = vector(mode = "list", length = numCuts)
      confl2 = vector(mode = "list", length = numCuts)
      pcs3 = vector(mode = "list", length = numCuts)
      confl3 = vector(mode = "list", length = numCuts)
      
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
        pcs2[[j]] = paste0(substr(partStrings[part], i, i+1), substr(partStrings[part], i+3, i+3))
        pcs3[[j]] = paste0(substr(partStrings[part], i, i), substr(partStrings[part], i+2, i+3))
        #So, this is tricky - it is the reverse complement of (1,3,4) that can conflict with (1,2,4) and vice versa
        #this needs to be taken care of later in the algorithm!
        confl2[[j]] = c(pcs2[[j]], complementary(reverse(pcs2[[j]])))
        confl3[[j]] = c(pcs3[[j]], complementary(reverse(pcs3[[j]])))
      }
      cutSeqsLoc[[part]] = cs
      piecesLoc[[part]] = pcs
      conflPiecesLoc[[part]] = confl
      piecesLoc2[[part]] = pcs2
      conflPiecesLoc2[[part]] = confl2
      piecesLoc3[[part]] = pcs3
      conflPiecesLoc3[[part]] = confl3
    }
    
    #so, I realize this is a bit odd, but practical in this case... (write to global vars)
    cutSeqs <<- cutSeqsLoc
    pieces <<- piecesLoc
    conflPieces <<- conflPiecesLoc
    pieces2 <<- piecesLoc2
    conflPieces2 <<- conflPiecesLoc2
    pieces3 <<- piecesLoc3
    conflPieces3 <<- conflPiecesLoc3
  } else {
    cutSeqs <<- NULL
    pieces <<- NULL
    conflPieces <<- NULL
    pieces2 <<- NULL
    conflPieces2 <<- NULL
    pieces3 <<- NULL
    conflPieces3 <<- NULL
  }
}

filterFindSeqCombData = function() {
  if (length(partStrings) != 0)
  {
    #remove cutSeqs that are not compatible with themselves (such as AATT) or that conflicts with the start and end tags
    startEndPieces = c(substr(startTag, 1, 3), substr(startTag, 2, 4), substr(endTag, 1, 3), substr(endTag, 2, 4) )
    startEndPieces2 = c(paste0(substr(startTag, 1, 2), substr(startTag, 4, 4)), paste0(substr(endTag, 1, 2), substr(endTag, 4, 4)))
    startEndPieces3 = c(paste0(substr(startTag, 1, 1), substr(startTag, 3, 4)), paste0(substr(endTag, 1, 1), substr(endTag, 3, 4)))

    cutSeqsLoc = cutSeqs
    piecesLoc = pieces
    conflPiecesLoc = conflPieces
    piecesLoc2 = pieces2
    conflPiecesLoc2 = conflPieces2
    piecesLoc3 = pieces3
    conflPiecesLoc3 = conflPieces3
    
    for (part in 1:length(partStrings)) {
      toRem = NULL
      numCuts = length(cutSeqsLoc[[part]])
      pcs = piecesLoc[[part]]
      confl = conflPiecesLoc[[part]]
      pcs2 = piecesLoc2[[part]]
      confl2 = conflPiecesLoc2[[part]]
      pcs3 = piecesLoc3[[part]]
      confl3 = conflPiecesLoc3[[part]]
      for (i in 1:numCuts) {
        #first (1,2,3) and (2,3,4) checks
        failure = any(pcs[[i]] %in% confl[[i]][3:4]) | any(startEndPieces %in% confl[[i]][1:4])#first check is compatible with self, second with start and end tags
        if (failure) {
          print(paste0("filtered 1 ", pcs[[i]][1], substr(pcs[[i]][2],3,3), " ",paste(pcs[[i]] %in% confl[[i]][3:4], collapse=" "), " - ",paste(startEndPieces %in% confl[[i]][1:4], collapse=" ")))
        }
        #then (1,2,4) need to compare (1,2,4) with rev compl of (1,3,4)
        failure = failure | pcs2[[i]] == confl3[[i]][2] | any(startEndPieces2 %in% confl2[[i]][1]) | any(startEndPieces2 %in% confl3[[i]][2])
        failure = failure | pcs3[[i]] == confl2[[i]][2] | any(startEndPieces3 %in% confl3[[i]][1]) | any(startEndPieces3 %in% confl2[[i]][2])
        if (failure) { 
          toRem = c(toRem, i)
          print(paste0("Removing: ", cutSeqsLoc[[part]][[i]]))
        } 
      }
      if (length(toRem) > 0) {
        print(-toRem)
        print(part)
        print(cutSeqsLoc)
        print(piecesLoc)
        print(conflPiecesLoc)
        cutSeqsLoc[[part]] = cutSeqsLoc[[part]][-toRem, drop=F]
        piecesLoc[[part]] = piecesLoc[[part]][-toRem, drop=F]
        conflPiecesLoc[[part]] = conflPiecesLoc[[part]][-toRem, drop=F]
      }
    }
    
    print("passed")
    
    #now write to global vars
    cutSeqs <<- cutSeqsLoc
    pieces <<- piecesLoc
    conflPieces <<- conflPiecesLoc
  }
}

checkData = function() {
  if (length(partStrings) != 0) {
    for (part in 1:length(partStrings)) {
      if (length(cutSeqs[[part]]) == 0) {
        return (FALSE)
      }
    }
  }  
  return (TRUE)
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
  confl2 = conflPieces2[[currPart]]
  confl3 = conflPieces3[[currPart]]
  numCuts = length(cutSeqs[[currPart]])
  for(i in 1:numCuts) {
    #check if this cut is compatible with the previous
    success = T
    if (!is.null(selVector)) {
      for (j in 1:length(selVector)) {
#        print (selVector)
#        print(j)
        success = success & !any(pieces[[j]][[selVector[j]]] %in% confl[[i]])
        #handle (1,2,4) - need to compare with rev compl of (1,3,4)
        success = success & (pieces2[[j]][[selVector[j]]] != confl2[[i]][1]) & (pieces2[[j]][[selVector[j]]] != confl3[[i]][2])
        #handle (1,3,4) - need to compare with rev compl of (1,2,4)
        success = success & (pieces3[[j]][[selVector[j]]] != confl3[[i]][1]) & (pieces3[[j]][[selVector[j]]] != confl2[[i]][2])
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

fwdFill = "gcat" #just space filler
revFill = "atgc" #just space filler
enz = "cgtctc" # recognition site for the restriction enzyme
#standardFwdEnd = "GTTTTAGAGCTAGAAATAGCAAGTTAA"
#standardRevEnd = "TTTCTTAGCTGCCTATACGGCAG"

genPrimers = function(positions, standardFwdEnd, standardRevEnd) {
  
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
      et = paste0(tolower(reverse(complementary(endTag))), "")
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
}


innerLoopTest = function(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) {
  success = TRUE
  for (j in 1:length(selVector)) {
    #        print (selVector)
    #        print(j)
    success = success & !any(pieces[[j]][[selVector[j]]] %in% confl[[i]])
    #handle (1,2,4) - need to compare with rev compl of (1,3,4)
    success = success & (pieces2[[j]][[selVector[j]]] != confl2[[i]][1]) & (pieces2[[j]][[selVector[j]]] != confl3[[i]][2])
    #handle (1,3,4) - need to compare with rev compl of (1,2,4)
    success = success & (pieces3[[j]][[selVector[j]]] != confl3[[i]][1]) & (pieces3[[j]][[selVector[j]]] != confl2[[i]][2])
  }
  
  return (success)
}

testStuff = function() {
  
  partStrings <<- c("ACGTTCAATTCCCGAAGCGG", "AAGACTTGTTAGCCC")
  numParts <<- length(partStrings)
  startTag <<- "CCAA"
  endTag <<- "GTTG"
  sideMargin <<- 3
  cutSeqs <<- NULL
  prepareFindSeqCombData()
  
  expRes = c("TTCC", "TCCC", "ATTC", "CCCG", "AATT", "CCGA", "CAAT", "CGAA", "TCAA", "GAAG", "TTCA")
  expRes2 = c("TTGT", "TGTT", "CTTG", "GTTA", "ACTT", "TTAG")
  #print(cutSeqs[[1]])
  #print(cutSeqs[[2]])
  #check cutSeqs
  print(all(cutSeqs[[1]] == expRes))
  print(all(cutSeqs[[2]] == expRes2))
  #check pieces (only check 1 item)
  all(pieces[[1]][[6]] == c("CCG", "CGA"))
  pieces2[[1]][[6]] == "CCA"
  pieces3[[1]][[6]] == "CGA"
  #check conflPieces (only check 1 item)
  all(conflPieces[[1]][[6]] == c("CCG", "CGA", "CGG", "TCG"))
  all(conflPieces2[[1]][[6]] == c("CCA", "TGG"))
  all(conflPieces3[[1]][[6]] == c("CGA", "TCG"))
  
  #now test filtering
  print(cutSeqs[[1]])# "TTCC" "TCCC" "ATTC" "CCCG" "AATT"(confl with itself) "CCGA"(confl start tag, (1,3,4)) "CAAT"(confl both start and end tag) 
  #"CGAA"(confl start tag(1,3,4)) "TCAA"(confl start tag) "GAAG" "TTCA"
  #start seq = "CCAA" (TTGG), end seq = "GTTG" (CAAC)
  print(cutSeqs[[2]])# "TTGT"(confl end tag) "TGTT"(confl end tag) "CTTG"(confl end tag) "GTTA"(confl end tag) "ACTT" "TTAG"(rev confl with start tag (1,2,4))
  filterFindSeqCombData()
  print(cutSeqs[[1]])#"TTCC" "TCCC" "ATTC" "CCCG" "GAAG" "TTCA"
  print(cutSeqs[[2]])#"ACTT"
  expRes = c("TTCC", "TCCC", "ATTC", "CCCG", "GAAG", "TTCA")
  expRes2 = c("ACTT")
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
  
  #now test the inner loop in findSeqComb, don't want to make it a function due to performance
  #So, the code is copied here instead. When testing, make sure the code is the same!

  #look at ATCG (CGAT)
  i=1
  confl= list(c("ATC","TCG","CGA","GAT"))
  confl2= list(c("ATG","CAT"))
  confl3= list(c("ACG","CGT"))
  
  #we also just check against 1 level of previous selections
  selVector = 1
  
  #TCGA should fail, first check
  pieces = list(list(c("TCG","CGA"))) #this is what we want to test against
  pieces2 = list(list(c("TCA"))) #this is what we want to test against
  pieces3 = list(list(c("TGA"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == FALSE

  #ATTG should fail, second check (1,2,4)
  pieces = list(list(c("ATT","TTG"))) #this is what we want to test against
  pieces2 = list(list(c("ATG"))) #this is what we want to test against
  pieces3 = list(list(c("ATG"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == FALSE
  
  #CGTT should fail, second check rev (1,2,4)
  pieces = list(list(c("CGT","GTT"))) #this is what we want to test against
  pieces2 = list(list(c("CGT"))) #this is what we want to test against
  pieces3 = list(list(c("CTT"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == FALSE
  
  #AACG should fail, third check (1,3,4)
  pieces = list(list(c("AAC","ACG"))) #this is what we want to test against
  pieces2 = list(list(c("AAG"))) #this is what we want to test against
  pieces3 = list(list(c("ACG"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == FALSE
  
  #CAAT should fail, third check rev (1,3,4)
  pieces = list(list(c("CAA","AAT"))) #this is what we want to test against
  pieces2 = list(list(c("CAT"))) #this is what we want to test against
  pieces3 = list(list(c("CAT"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == FALSE

  #Now test some things that should work
  #CAGT should work
  pieces = list(list(c("CAG","AGT"))) #this is what we want to test against
  pieces2 = list(list(c("CAT"))) #this is what we want to test against
  pieces3 = list(list(c("CGT"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == TRUE
  
  #GGAA should work
  pieces = list(list(c("GGA","GAA"))) #this is what we want to test against
  pieces2 = list(list(c("GGA"))) #this is what we want to test against
  pieces3 = list(list(c("GAA"))) #this is what we want to test against
  
  innerLoopTest(selVector, pieces, pieces2, pieces3, confl, confl2, confl3, i) == TRUE
  
  
}



onGenButton = function(input, output, session) {
  strWarn = ""
  res = setGlobalVars(input)
  if (res[[1]]) {
    prepareFindSeqCombData()
    filterFindSeqCombData()
    if (checkData()) {
    
      # All global variables set, run the recursive function to find the place to cut each part
      res = findSeqComb()
      if (res[[1]]) {
      
        #extract the resulting sticky end sequences
        seqs = extractSeqs(res)
    
        #now find where this is in the parts:
        positions = extractPositions(seqs)
        print(positions)
        
        #generate primers:
        primers = genPrimers(positions, input$FwdBP, input$RevBP)
        dfPrim = data.frame(Pair=1:length(primers[[1]]), Forward = primers[[1]], Reverse = primers[[2]], stringsAsFactors = F)
        colnames(dfPrim) = c("Pair", "Forward primer", "Reverse primer")
        
        output$outpPrim = DT::renderDataTable(DT::datatable({dfPrim}, escape = F, width="100%", rownames = FALSE, class="compact", selection = 'none', options = list(dom = 't', pageLength = 20)))
        
        allSeqs = c(startPart, partStrings, endPart)
          partStrings = NULL
        startPart = NULL
        endPart = NULL
        
        for (s in 1:length(allSeqs)) {
          for (i in 1:length(enzymeRecSitesIds)) {
            ind = strfind(allSeqs[[s]], enzymeRecSitesSeqs[i])
            if (length(ind) > 0) {
              strWarn = paste0(strWarn, "Restr. enzyme recogn. site found in seq: ", enzymeRecSitesIds[i], ": ", enzymeRecSitesSeqs[i], " in part ", s, "\n")
            }
            #if the enzymes are not the same after reverse and complement (which some are), look for the reverse complement as well
            #This was tested manually
            revCompl = reverse(complementary(enzymeRecSitesSeqs[i]))
            if (enzymeRecSitesSeqs[i] != revCompl) {
              ind = strfind(allSeqs[[s]], revCompl)
              if (length(ind) > 0) {
                strWarn = paste0(strWarn, "Restr. enzyme recogn. site found in the reverse compl. seq: ", enzymeRecSitesIds[i], ": ", revCompl, " in part ", s, "\n")
              }
            }
          }
        }
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
  if (strWarn == "") {
    output$warningsTitle = renderText("")
  } else {
    output$warningsTitle = renderText({HTML("<b>Warnings</b>")})
  }
  output$warnings = renderText(strWarn)
}

server <- function(input, output, session) {
  observeEvent(input$genButton, {
    onGenButton(input, output, session)
  })
  
}



#don't add any code after the shinyApp line below, it has to be the last line!
shinyApp(ui = ui, server = server)
