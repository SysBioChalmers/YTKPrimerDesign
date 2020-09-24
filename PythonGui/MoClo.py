import tkinter as tk
from tkinter import ttk
import tkinter.font as tkFont
import primer3 as p3
import tkinter.filedialog as tkFile
import os


#print(p3.calcTm('GTAAAACGACGGCCAGT'))



def reverse(s):
    return s[::-1]

def compChar(s):
    compSwitcher = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }
    return compSwitcher.get(s, "X")

print(compChar("A"))

def cleanupSeq(s):
    s = s.upper()
    ret = ""
    for i in range(len(s)):
        if (s[i] == 'A') | (s[i] == 'T') | (s[i] == 'C') | (s[i] == 'G'):
            ret += s[i]
    return ret

def complementary(s):
    ret = ""
    for i in range(len(s)):
        ret += compChar(s[i])
    return ret

def gcContent(s):
    gcs = 0
    for i in range(len(s)):
        if (s[i] == 'G') | (s[i] == 'C'):
            gcs += 1
    return gcs/len(s)

#test
print(gcContent("AATTGGCCAATTGGCC"))


def meltTemp(s):
    return (p3.calcTm(s))
def hairpin(s):
    return (p3.calcHairpin(s))

#print(hairpin('AATTTCCCCGGGG'))

def checkPrimer(s):
    gcCont = gcContent(s)
    #todo: add borders for gcContent?
    mt = meltTemp(s)
    if (mt < 48) | (mt > 60):
        return False
    return True
print(complementary("ATTCG"))
print(gcContent("ATTCG"))

#for forward primer
univStartTag = "GCATCGTCTCATCGGTCTCA"
startTags = dict()
#startTags["Type 1"] = "CCCT"
startTags["Type 2"] = "AACG"
startTags["Type 3"] = "TATG"
startTags["Type 3A"] = "TATG"
startTags["Type 3B"] = "TTCT"
startTags["Type 4"] = "ATCCTAACTCGAG"
startTags["Type 4A"] = "ATCC"
startTags["Type 4B"] = "TGGC"
#startTags["Type 5"] = "GCTG"
startTags["Type 6"] = "TACA"
startTags["Type 7"] = "GAGT"
startTags["Type 8"] = "CCGAGCGGCCGC"
startTags["Type 8A"] = "CCGAGCGGCCGC"
startTags["Type 8B"] = "CAAT"

univEndTag = "TGAGACCTGAGACGGCAT"

endTags = dict()
#endTags["Type 1"] = "AACG"
endTags["Type 2"] = "AGATCTATG"
endTags["Type 3"] = "GGATCC"
endTags["Type 3A"] = "GGTTCT"
endTags["Type 3B"] = "GGATCC"
endTags["Type 4"] = "GCTG"
endTags["Type 4A"] = "TAACTCGAGTGGC"
endTags["Type 4B"] = "GCTG"
#endTags["Type 5"] = "TACA"
endTags["Type 6"] = "GAGT"
endTags["Type 7"] = "CCGA"
endTags["Type 8"] = "GCGGCCGCCCCT"
endTags["Type 8A"] = "GCGGCCGCCAAT"
endTags["Type 8B"] = "CCCT"


enzymeRecSitesIds = ["BsaI", "BsmBI", "NotI", "BglII", "BamHI", "XhoI"]
enzymeRecSitesSeqs = ["GGTCTC", "CGTCTC", "GCGGCCGC", "AGATCT", "GGATCC", "CTCGAG"]



root = tk.Tk()
# s = ttk.Style()
# s.theme_use('clam')

root.title('MoClo Primer Design Tool')
root.minsize(400, 300)
inp = tk.Text(root, width=140, height=10)
inp.focus()

outpFwd = tk.Text(root, width=140, height=12)
outpFwd.config(state='disabled')

outpRev = tk.Text(root, width=140, height=12)
outpRev.config(state='disabled')

#make font a bit larger
fontStyle = tkFont.Font(family="Arial", size=14)

filterVar = tk.IntVar(value=1)
filterChk = tk.Checkbutton(root, text="Filter primers", variable=filterVar, height=3, font=fontStyle)

warnings = tk.Text(root, width=70, height=8, foreground="red")
warnings.config(state='disabled')


def onGenButton():
    #expButton.config(state='normal') not needed
    seq = cleanupSeq(inp.get("1.0", tk.END))
    # loop through all primer lengths
#    forwardPrimers = list()
#    reversePrimers = list()
    startId = startCombo.get()
    startTag = startTags[startId]
    endId = endCombo.get()
    endTag = endTags[endId]
    fwInd = 1
    revInd = 1

    outpFwd.config(state='normal')
    outpFwd.delete("1.0", tk.END)
    outpRev.config(state='normal')
    outpRev.delete("1.0", tk.END)

    global fwdPrimers
    global revPrimers
    global fwdPrimersText
    global revPrimersText
    fwdPrimers = list()
    revPrimers = list()
    fwdPrimersText = list()
    revPrimersText = list()



    for length in range(17, 29):
        forwardPrimer = seq[0:length]
        if (filterVar.get() == 0) | checkPrimer(forwardPrimer):
            fwdPrimers.append(fwdPrimers)
            outputString = univStartTag + startTag + forwardPrimer + '\tGC cont: ' + '{:d}'.format(round(gcContent(forwardPrimer)*100))+ '%' + '\tLength: ' + str(len(forwardPrimer)) + '\tMelt temp: ' + '{:.01f}'.format(meltTemp(forwardPrimer)) + ' Hairpin Tm: {:.01f}'.format(hairpin(forwardPrimer).tm) + '\n'
            fwdPrimersText.append(outputString)
            outpFwd.insert(str(fwInd) + ".0", outputString)
            fwInd += 1

        reversePrimer = reverse(seq[-length::]) #reverses automatically
        if (filterVar.get() == 0) | checkPrimer(reversePrimer):
            revPrimers.append(reversePrimer)
            outputString = complementary(reverse(univEndTag)) + endTag + complementary(reversePrimer) + '\tGC cont: ' + '{:d}'.format(round(gcContent(reversePrimer)*100))+ '%' + '\tLength: ' + str(len(reversePrimer)) + '\tMelt temp: ' + '{:.01f}'.format(meltTemp(reversePrimer)) + ' Hairpin Tm: {:.01f}'.format(hairpin(reversePrimer).tm)  + '\n'
            revPrimersText.append(outputString)
            outpRev.insert(str(revInd) + ".0", outputString)
            revInd += 1
    #remove the last newline
    outpFwd.delete(str(fwInd) + ".0", tk.END)
    outpRev.delete(str(revInd) + ".0", tk.END)
    outpFwd.config(state='disabled')
    outpRev.config(state='disabled')

    #add warnings:
    #enzyme recognition sites:
    warnings.config(state='normal')
    warnings.delete("1.0", tk.END)

    warnInd = 1

    for i in range(len(enzymeRecSitesIds)):
        a = -1
        while True:
            a = seq.find(enzymeRecSitesSeqs[i], a+1)
            if a == -1:
                break
            warnings.insert(str(warnInd) + ".0", "Restr. enzyme recogn. site found in seq: " + enzymeRecSitesIds[i] + ": " + enzymeRecSitesSeqs[i] + " at pos " + str(a) + "\n")
            warnInd += 1
    #issue warnings for certain parts, since the start/end codon is included in the primer
    if (startId == "Type 3") | (startId == "Type 3A"):
        warnings.insert(str(warnInd) + ".0",
                        "Start codon is already included in fwd primer. Do not include stop codon in pasted sequence if part will be used as fusion or tagged protein.\n")
        warnInd += 1
    if endId == "Type 3B":
        warnings.insert(str(warnInd) + ".0",
                        "Stop codon is already included in rev primer. Do not include stop codon in pasted sequence if part will be used as fusion or tagged protein.\n")
        warnInd += 1

    warnings.delete(str(warnInd) + ".0", tk.END)
    warnings.config(state='disabled')

#win = None #just assign to something


def onExpButton():
    #global win
    #win = tk.Toplevel()

    #win.wm_title("Window")
    #win.grab_set()
    #lb = tk.Listbox(win)
    #lb.pack
    #but = tk.Button(win, command=onFinishExp)
    #but.pack();


    filename = tkFile.asksaveasfilename(initialdir=os.getcwd(), title="Select file", defaultextension='.FASTA', filetypes=(("fasta files","*.FASTA"),("all files","*.*")))
    print(filename)

    fwdOverhangStd = "TCGGTCTCA"
    fwdOverhangTypeSpec = startTags[startCombo.get()]
    seq = cleanupSeq(inp.get("1.0", tk.END))
    revOverhangStd = "TGAGACC"
    revOverhangTypeSpec = endTags[endCombo.get()]
    backbone = "AGACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCAGATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCCAAGCGAGCTCGATATCAAATTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAAACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGAAATTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCCGTAATATCCAGCTGAACGGTCTGGTTATAGGTACATTGAGCAACTGACTGAAATGCCTCAAAATGTTCTTTACGATGCCATTGGGATATATCAACGGTGGTATATCCAGTGATTTTTTTCTCCATTTTAGCTTCCTTAGCTCCTGAAAATCTCGATAACTCAAAAAATACGCCCGGTAGTGATCTTATTTCATTATGGTGAAAGTTGGAACCTCTTACGTGCCCGATCAATCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTAG"

    totSeqToExp = fwdOverhangStd + fwdOverhangTypeSpec + seq + revOverhangTypeSpec + revOverhangStd + backbone

    #now write the FASTA file:
    file1 = open(filename, "w")
    file1.write(">Yeast toolkit export. Starts at " + startCombo.get() + ", ends at " + endCombo.get() + "\n")
    ll = len(totSeqToExp)
    for i in range(0, ll, 80):
        file1.write(totSeqToExp[i:min(i+80, ll)] + '\n')

    file1.close()


#def onFinishExp():
#    win.destroy()


buttonFrame = tk.Frame(root)

genButton = tk.Button(buttonFrame,
                   command=onGenButton,
                   text="Generate\nprimers",
                   width=13,
                   height=2,
                   font=fontStyle
)

expButton = tk.Button(buttonFrame,
                   command=onExpButton,
                   text="Export\nFASTA",
                   width=13,
                   height=2,
                   font=fontStyle
)
#expButton.config(state='disabled') #not needed

#test = tk.Listbox(root)
#test.insert(1, "Nachos")
#test.insert(2, "Sandwich")
#test.insert(3, "Burger")
#test.insert(4, "Pizza")
#test.insert(5, "Burrito")


startCombo = ttk.Combobox(root, values=list(startTags.keys()), font=fontStyle)

#comboExample.grid(column=0, row=1)
startCombo.config(state='readonly')
#startCombo.pack()
startCombo.current(0)

endCombo = ttk.Combobox(root, values=list(endTags.keys()), font=fontStyle)

#comboExample.grid(column=0, row=1)
endCombo.config(state='readonly')
#endCombo.pack()
endCombo.current(0)

root.option_add("*TCombobox*Listbox*Font", fontStyle)



#set up a grid layout
tk.Label(root, text="Input sequence:", font=fontStyle).grid(row=0, column=0)
inp.grid(row=0, column=1, columnspan=4)

tk.Label(root, text="Part type prefix:", height=2, font=fontStyle).grid(row=1, column = 1)
startCombo.grid(row=1, column=2, sticky=tk.W+tk.E)
tk.Label(root, text="Part type suffix:", height=2, font=fontStyle).grid(row=2, column=1)
endCombo.grid(row=2, column=2, sticky=tk.W+tk.E)
filterChk.grid(row=3, column=1)
buttonFrame.grid(row=3, column=2, sticky=tk.W+tk.E)

#test.grid(row=1, column=0, columnspan=3)

tk.Label(root, text="Warnings:", font=fontStyle).grid(row=1, column=3, rowspan=3, sticky=tk.N+tk.S)
warnings.grid(row=1, column=4, rowspan=3, sticky=tk.W+tk.E+tk.N+tk.S)

tk.Label(root, text="Forward primers:", font=fontStyle).grid(row=4, column=0)
outpFwd.grid(row=4, column=1, columnspan=4)
tk.Label(root, text="Reverse primers:", font=fontStyle).grid(row=5, column=0)
outpRev.grid(row=5, column=1, columnspan=4)

#the button frame layout
genButton.grid(row=0, column=0, sticky=tk.W+tk.E)
expButton.grid(row=0, column=1, sticky=tk.W+tk.E)


root.mainloop()


