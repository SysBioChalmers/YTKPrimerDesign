

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




partStrings = c("CAACAAGTTATTACTCAAA",
          "TTGGTTGATGATAGATTCA",
          "TTCAAATCTCTGATTCTAT",
          "AGAAATCTGAAGGTTTGGC",
          "CTACTGATTCTACAAAAAG",
          "ATCACAAGAACAACCAATC",
          "GATAAGGATCCAATTAAAG",
          "GCAACTGCTGCAATGGCAA",
          "ACACCATTAGTTAAGGAAC",
          "CAAGATACATGGTATTACC",
          "CACCAGATATCGCTAACGA",
          "TGCAATCAATTAATTTGCC")


#parts = vector(mode = "list", length = length(partStrings))
#for (part in 1:length(partStrings)) {
#  parts[[part]] = strsplit(partStrings[part],"")[[1]]
#}

#could get conflict either with other sequences or reverse complement of other sequences (including self, should be checked at once)
#enough to have an overlap of 3 nucleotides to create a conflict (unclear exactly how)

#first create all possible sequences and their conflicting pieces of length 3
cutSeqs = vector(mode = "list", length = length(partStrings)) #the 4-letter cut site
pieces = vector(mode = "list", length = length(partStrings)) #the 2 three letter pieces
conflPieces = vector(mode = "list", length = length(partStrings)) #the 4 three letter pieces that needs to be checked against the pieces of other cutSeqs
for (part in 1:length(parts)) {
  numCuts = str_length(partStrings[part])-3
  cs = vector(mode = "list", length = numCuts)
  pcs = vector(mode = "list", length = numCuts)
  confl = vector(mode = "list", length = numCuts)
  for (i in 1:numCuts) {
    cs[[i]] = substr(partStrings[part], i, i+3)
    pcs[[i]] = c(substr(partStrings[part], i, i+2), substr(partStrings[part], i+1, i+3))
    confl[[i]] = c(substr(partStrings[part], i, i+2), substr(partStrings[part], i+1, i+3), 
                       complementary(reverse(substr(partStrings[part], i, i+2))), complementary(reverse(substr(partStrings[part], i+1, i+3))))
  }
  cutSeqs[[part]] = cs
  pieces[[part]] = pcs
  conflPieces[[part]] = confl
}

#remove cutSeqs that are not compatible with themselves (such as AATT)
for (part in 1:length(parts)) {
  toRem = NULL
  numCuts = length(cutSeqs[[part]])
  pcs = pieces[[part]]
  confl = conflPieces[[part]]
  for (i in 1:numCuts) {
    if (any(pcs[[i]] %in% confl[[i]][3:4])) {
      toRem = c(toRem, i)
      print(paste0("Removing: ", cutSeqs[[part]][[i]]))
    }
  }
  if (length(toRem) > 0) {
    cutSeqs[[part]] = cutSeqs[[part]][-toRem]
    pieces[[part]] = pieces[[part]][-toRem]
    conflPieces[[part]] = conflPieces[[part]][-toRem]
  }
}

#TODO: Remove cutSeqs that are not compatible with other sticky ends present

#strategy:
#go through the parts in order
#loop through all cutSeqs for each part and check if that works with the selected cutSeqs in the previous parts
#for each cutSeq that works, select that and try to proceed recursively with the next part. If it fails, continue the loop.

#returns a list, where the first value in the list is true if it succeeded, and the second is the index vector
numParts = length(parts)
findSeqComb = function(selVector) {
  currPart = length(selVector) + 1
#  if (currPart == 6) {
#    print("test")
#  }
  if (currPart > numParts) {
    return (list(T, selVector)) #we have succeeded, there are no more to find
  }
  confl = conflPieces[[currPart]]
  numCuts = length(cutSeqs[[part]])
  for(i in 1:numCuts) {
    #check if this cut is compatible with the previous
    success = T
    for (j in selVector) {
      success = success & !any(pieces[[j]][[selVector[j]]] %in% confl[[i]])
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

findSeqComb(NULL)
  
  
  
  
