# -*- coding: utf-8 -*-

# To generate permutations
from itertools import product 

#For GUI
from tkinter import Tk, mainloop, LEFT, TOP , scrolledtext
from tkinter.ttk import *
from tkinter import *
import tkinter
import tkinter.font as tkFont
import string




#Utility function to convert tuple to string    
def convertTuple(tup): 
    str =  ''.join(tup) 
    return str

'''
Class is defined for Instruction text formatting 
'''
class RichText(Text):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        default_font = tkFont.nametofont(self.cget("font"))

        em = default_font.measure("m")
        default_size = default_font.cget("size")
        bold_font = tkFont.Font(**default_font.configure())
        italic_font = tkFont.Font(**default_font.configure())
        h1_font = tkFont.Font(**default_font.configure())
        normal_font = tkFont.Font(**default_font.configure())
        normal_u_font = tkFont.Font(**default_font.configure())
        b_u_font = tkFont.Font(**default_font.configure())
        
        bold_font.configure(size= 16,family="Arial",weight="bold")
        italic_font.configure(size= 16,family="Arial",slant="italic")
        h1_font.configure(size= 18,family="Arial", weight="bold",underline =-1)
        normal_font.configure(size= 16,family="Arial")
        normal_u_font.configure(size= 16,family="Arial",underline =-1)
        b_u_font.configure(size= 16,family="Arial", weight="bold" ,underline =-1)

        self.tag_configure("bold", font=bold_font)
        self.tag_configure("italic", font=italic_font)
        self.tag_configure("noraml", font=normal_font)
        self.tag_configure("h1", font=h1_font, spacing3=default_size)
        self.tag_configure("boldunderline", font = b_u_font)
        self.tag_configure("normalunderline", font = normal_u_font)
        
 

        lmargin2 = em + default_font.measure("\u2022 ")
        self.tag_configure("bullet", lmargin1=em, lmargin2=lmargin2)

    def insert_bullet(self, index, text):
        self.insert(index, f"\u2022 {text}", "bullet")

class REMP:
    
    def __init__(self, table1,table2):
        self.table1 = table1
        self.table2 = table2
        
                
    def generateDegenerateSequences(self,nu_seq):
        
        # Generating triplets from Nucleotide sequence
        nu_seq = nu_seq.translate(str.maketrans('', '', ' \n\t\r')) #Removing White spaces
        nu_seq = nu_seq.upper()
        
        codons_list = {}        
        no_of_codons = int(len(nu_seq)/3)
        nu_seq_trim = nu_seq[:no_of_codons*3]
        print("Input Sequence \n",nu_seq_trim)
        # print("no_of_codons", no_of_codons)
        i = 0
        while(i < no_of_codons*3):
            codons_list['condon_'+str(int(i/3))] = [nu_seq[i:i+3]]
            i += 3
        
        # Identifying substitute triplets for given Nucleotide sequence
        for codon in codons_list:
            for triplet in self.table1:
                if codons_list[codon][0] in self.table1[triplet]:
                    for sub in self.table1[triplet] :
                        if sub not in codons_list[codon]:
                            codons_list[codon].append(sub)                            
                  
        #Generating substitute triplets Sequences after substitute triplets       
        codons_substitut_lists = []
        for codon in codons_list:
            codons_substitut_lists.append(codons_list[codon])
        # print("codons_substitut_list",codons_substitut_lists)
        
        # Window contains 4 triplets, restriction site can obtain from 4 triplets 
        window_size = 4
        no_windows = 0
        if no_of_codons < 4:
            no_windows = 1
        else:
            no_windows = no_of_codons - (window_size-1)
        
        # print("no_windows ",no_windows)
               

        
        '''
        Begin Sliding window technique
        
        Window contains 4 triplets , next iteration from previous window right most triplet removed ,
        next triplet after left most of previous window will be added in current window. Will generate the
        all possible degenerate sequence within  window.
        
        Prefix and Postfix of sequence kept same.
        
        '''
        result= {}
        match_count = 0
        for wind in range(0,no_windows):
            window_conodns_list = codons_substitut_lists[wind:wind+4]
            
            # Generating  degenerate sequences within window
            windo_degenerate_sequences = list(product(*window_conodns_list))
            windo_degenerate_sequences_list = [] 
            for tup in windo_degenerate_sequences:
                windo_degenerate_sequences_list.append(convertTuple(tup))
            
            # Identifying the restriction sites in degenerate sequences, each restriction site can obtain from multiple sequences. 
            for Palindrome in self.table2:
                    for seq in windo_degenerate_sequences_list:
                        ind = seq.find(Palindrome)
                        if ind != -1:
                            match_count += 1
                            prefix_ind =  wind*3
                            posfix_ind = prefix_ind+12
                            # Adding prefix and postfix of given sequencre to degenarate sequnces 
                            tot_seq = nu_seq_trim[:prefix_ind]+seq+nu_seq_trim[posfix_ind:]
                            if Palindrome not in result:
                                result[Palindrome] = [tot_seq]                               
                            else:
                                result[Palindrome].append(tot_seq)

        '''
        End Sliding window technique

        '''
        
        #Identifying minimum no of alphabets changed degenerate sequence for each restriction site
        uniq_min_seq = []
        uniq_enzymes = []
        for hexamer in result:
            if hexamer not in uniq_enzymes:
                uniq_enzymes.append(hexamer)
            min_changes = len(nu_seq_trim)
            min_seq = ""
            for degen_seq in result[hexamer]:
                changes_count = 0
                for ind in range(0,len(degen_seq)):
                    if (degen_seq[ind] != nu_seq_trim[ind]):
                        changes_count += 1
                if changes_count <  min_changes:
                    min_changes = changes_count
                    min_seq = degen_seq
                if changes_count == len(hexamer):
                    ind = len(degen_seq)
            if min_seq not in uniq_min_seq:
                uniq_min_seq.append(min_seq)
        
        # Sorting the all degenerate sequences based on minimum no of alphabets changed . 
        for i in range(len(uniq_min_seq)):
            min_ind = -1
            min_changes = len(uniq_min_seq[i])+1
            for j in range(i,len(uniq_min_seq)):
                changes_count = 0
                for ind in range(0,len(uniq_min_seq[j])):
                    if (uniq_min_seq[j][ind] != nu_seq_trim[ind]):
                        changes_count += 1
                if changes_count < min_changes:
                    min_changes = changes_count
                    min_ind = j
            uniq_min_seq[i],uniq_min_seq[min_ind] = uniq_min_seq[min_ind],uniq_min_seq[i]
        
        uniq_seq_enzymes = {}
        
        # Final output for the given input sequence
        for uniq_seq in uniq_min_seq:
            for uniq_enz in uniq_enzymes:
                if uniq_seq.find(uniq_enz) != -1:
                    if uniq_seq not in uniq_seq_enzymes:
                        uniq_seq_enzymes[uniq_seq] = [ ' '.join(x+' ('+uniq_enz+')' for x in self.table2[uniq_enz])]
                    else:
                        uniq_seq_enzymes[uniq_seq].append(' '.join(x+' ('+uniq_enz+')' for x in self.table2[uniq_enz]))
        print("Final Output Sequences with Enzymes\n",uniq_seq_enzymes)
        
        return uniq_seq_enzymes
              
       
if __name__=='__main__':

    #Loading Degenerate triplets from file table1.txt 
    table1 = {}    
    with open('./table1.txt', 'r') as document:
        for line in document:
            if line.strip():  
                key, value = line.split(None, 1)
                value = value.strip()
                table1[key] = value.split(',')
    
    #Loading Restriction site and its enzyme from file table2.txt
    table2 = {}
    with open('./table2.txt', 'r') as document:        
        for line in document:
            if line.strip():  
                key, value = line.split(None, 1)  
                value = value.strip()
                table2[key] = value.split(',')
                    
    # print("Table 1 \n",table1)
    # print("Table 2 \n",table2)
          
    #UI Part
    '''
    Instructions regarding  input and output of the program.
    '''
    def clickInstructions():
        
        instru_window = Toplevel(root)
        instru_window.title("Instructions")
        instru_window.geometry("1100x650") 
        bg_color = "ivory3"
        instru_window.configure(background=bg_color)
        
        Tops = Frame(instru_window, width = 800, relief = FLAT,bg = bg_color)
        Tops.pack(side = TOP,expand=True, fill=BOTH)
                
        text = RichText(Tops, width=50, height=15,padx = 70)
        yscrollbar = Scrollbar(Tops, orient=VERTICAL, command=text.yview)
        text['yscroll'] = yscrollbar.set
        yscrollbar.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        
        text.insert("end", "                                                   ","noraml")
        text.insert("end", "Instructions", "h1")
        text.insert("end", "                             \n\n","noraml")
    
        text.insert("end", "Steps involved in using REMP software are explained below using an example:\n", "bold")
        text.insert("end", "To mutate amino acid Glutamate (E) in the following ORF to a Lysine (K).\n","noraml")
        text.insert("end", "GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "GAA ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC.\n","noraml")
        text.insert("end", "  A       D      S      D     P     L    ","noraml")
        text.insert("end", " E    ","bold")
        text.insert("end", "  S       A      D      V     S     D  \n","noraml")
        text.insert("end","STEP 1","boldunderline")
        text.insert("end"," Generate a mutant sequence by replacing GAA codon of E with ", "noraml")
        text.insert("end","any codon ", "bold")
        text.insert("end", "of Lysine.\n","noraml")
        text.insert("end", " GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAA ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC [OR]\n","noraml")
        text.insert("end", " GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC \n","noraml")
        text.insert("end", "   A       D      S      D     P     L    ","noraml")
        text.insert("end", " K    ","bold")
        text.insert("end", "   S       A      D      V     S     D  \n","noraml")
        text.insert("end","____________________________________________________________________________________________________________________\n\n","normal")
        
        text.insert("end","STEP 2","boldunderline")
        text.insert("end", " Input the mutant sequence to REMP software. In the input sequence, include ","noraml")
        text.insert("end", "complete \n","bold")
        text.insert("end", "codon ","bold")
        text.insert("end", "of the amino acid at 5` end and codons of ","noraml")
        text.insert("end", "at least 4 amino acids on either side ","normalunderline")
        text.insert("end", "of intended \n","noraml")
        text.insert("end", "mutation (for efficient amplification in PCR).\n\n","noraml")
        
        text.insert("end", "TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG [4 codons on sides]\n","noraml")
        text.insert("end", "GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA [5 codons on sides]\n","noraml")
        text.insert("end", "GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC [6 codons on sides]\n\n","noraml")
        text.insert("end","Wrong ","bold")
        text.insert("end", "input sequences do not have complete codon of amino acid at 5` end. The codon of amino\n","noraml")
        text.insert("end", "acid at 3` end is not required to be a complete one.\n","noraml")   
        text.insert("end","____________________________________________________________________________________________________________________\n\n","normal")
        
        text.insert("end","STEP 3","boldunderline")
        text.insert("end", " From the output sequences, generate a mutant primer with a desired restriction site for\n","noraml")
        text.insert("end", "screening the mutation. An ","noraml")
        text.insert("end", "ideal mutant primer ","bold")
        text.insert("end", "is one with silent mutations as ","noraml")
        text.insert("end", "close to intended\n","bold")
        text.insert("end", "mutation ","bold")
        text.insert("end","(E to K in above example) as possible, and, with a ","noraml")
        text.insert("end","least number of base changes ","bold")
        text.insert("end","from\n","noraml")
        text.insert("end","the input sequence, if possible.\n\n","noraml")
        text.insert("end","    - ","bold")
        text.insert("end","REMP output displays mutant sequences carrying a restriction site(s) created via \n","noraml")
        text.insert("end","      silent mutation(s).\n","bold")
        text.insert("end","    - ","bold")
        text.insert("end","The output sequences are sorted according to number of changes from input sequence, with\n","noraml")
        text.insert("end","      the least changes being on ","noraml")
        text.insert("end","top of the list ","bold")
        text.insert("end","and the most at the bottom of the list.\n","noraml")
        text.insert("end","    - ","bold")
        text.insert("end","The name of the restriction enzyme along with its recognition site is displayed at the right\n","noraml")
        text.insert("end","      end of the output sequence.\n\n","noraml")
        text.insert("end","TCAGACCC","noraml")
        text.insert("end","T","boldunderline")
        text.insert("end","TTAAAGAGCGCCGACGTG\t\t\t\t\t\t\t","noraml")
        text.insert("end","DraI (TTTAAA)\n","noraml")
        text.insert("end","TCAGACCCATTAAAGAGCGC","noraml")
        text.insert("end","T","boldunderline")
        text.insert("end","GACGTG\t\t\t\t\t\t\t","noraml")
        text.insert("end","AfeI (AGCGCT)\n\n","noraml")
        text.insert("end","Of the two outputs shown above, both DraI and AfeI are generated via 1 base change, however,\n","noraml")
        text.insert("end","DraI is closest to intended mutation of E(GAA) to K(AAG) and would be ideal for a primer design.\n","noraml")
        text.insert("end","____________________________________________________________________________________________________________________\n\n","normal")
        text.insert("end","NOTES:","boldunderline")
        text.insert("end"," (1) The output sequence can be copied and final primer with required Tm can be designed\n","noraml")
        text.insert("end"," using a software such as ‘Oligo Calc’ ","noraml")
        text.insert("end","[http://biotools.nubic.northwestern.edu/OligoCalc.html].\n\n","boldunderline")
        text.insert("end","(2) A restriction site generating distinguishable digested DNA pattern of mutant plasmid compared\n","noraml")
        text.insert("end","to a wild type plasmid, can be analyzed using a software such as NEBcutter V2.0\n","noraml")
        text.insert("end","[http://nc2.neb.com/NEBcutter2/].\n","boldunderline")
        text.configure(state ='disabled')

    # Functionality based on input from UI buttons.       
    def doREMP():
        # Getting input sequence from textbox 
        nu_seq = input_text.get().strip()
        # print("Nucleotide sequence= ", (input_text.get()))
        outputText.configure(state ='normal')
        outputText.delete(1.0, END)
        rempObj = REMP(table1,table2)
        uniq_seq_enzymes = rempObj.generateDegenerateSequences(nu_seq)
        xscrollbar.config(command=outputText.xview)
        yscrollbar.config(command=outputText.yview)
        outputText.configure(state ='normal')
        result = ""
        if (len(uniq_seq_enzymes) == 0):
            outputText.insert('end -1 chars',"None")
        else:
            for seq in uniq_seq_enzymes:
                result += seq+'\t\t\t\t\t'+' '.join(x for x in uniq_seq_enzymes[seq])
                result += '\n\n'
                                
            outputText.insert('end -1 chars',result)
        outputText.configure(state ='disabled')
    
    # Clear the input and outputs
    def Clear():        
        input_text.set("")
        outputText.configure(state ='normal')
        outputText.delete(1.0, END)
        outputText.configure(state ='disabled')
    
    
    # Tkinter GUI
    root = Tk()
    root.title("REMP")
    root.geometry("1200x600")
    bg_color = "gray"
    root.configure(background=bg_color)
    

    Tops = Frame(root, width = 1100, relief = FLAT,bg=bg_color)
    Tops.pack(side = TOP) 
              
    lblTitle = Label(Tops, font = ('Arial', 40, 'bold'), 
                    text = "REMP ", fg = "cyan2", bd = 10, bg=bg_color, justify = LEFT)
    lblTitle.grid(row = 0, column = 0)
    
    lblTitle2 = Label(Tops, font = ('Arial', 20, 'bold'), 
                    text = "REstriction site in Mutant Primer", fg = "cyan2", bd = 10, bg=bg_color, justify = LEFT)
    lblTitle2.grid(row = 0, column = 1,columnspan = 3)    
    
    lblInput = Label(Tops, font = ('Arial', 18), bg=bg_color,
                text = "Mutant Primer",fg = "thistle1", bd = 5, justify=RIGHT)
    lblInput.grid(row = 1, column = 0)
    
    lblInput2_txt = '[Input complete codon of amino acid at 5` end]'
    lblInput = Label(Tops, font = ('Arial', 16), bg=bg_color,
                text = lblInput2_txt,fg = "thistle1", bd = 5, justify=LEFT)
    lblInput.grid(row = 1, column = 1)
    
    input_text = StringVar()

    txtInput = Entry(Tops,  font = ("Arial",16),
                         bd = 2,  width = 90,
                         textvariable = input_text)                          
    txtInput.grid(row = 2, column = 0,columnspan = 5,ipady=4)
            

    btnInstru = Button(Tops,  padx = 2, pady = 2, bd = 2, fg = "yellow", 
                        font = ('Arial', 16), width = 12, 
                        text = "INSTRUCTIONS",relief = FLAT,bg = bg_color,
                command=clickInstructions,justify = LEFT).grid(row = 1, column = 4)
    
    
    btnSubmit = Button(Tops, padx = 2, pady = 2, bd = 2, 
                  fg = "green2", font = ('Arial', 16), 
                    width = 8, text = "SUBMIT", bg = bg_color,relief = FLAT,justify = RIGHT,
                    command = doREMP).grid(row = 3, column = 0)

    btnClear = Button(Tops, padx = 2, pady = 2, bd = 2, 
              fg = "OliveDrab1", font = ('Arial', 16), 
                width = 8, text = "CLEAR", bg = bg_color,relief = FLAT,justify = RIGHT,
                command = Clear).grid(row = 3, column = 4)
    

    
    lblOutput = Label(Tops, font = ('Arial', 18), bg=bg_color,
                text = "Mutant Primer With a Restriction Site", bd = 5,fg = "thistle1",justify=LEFT)
    lblOutput.grid(row = 4, column = 0,columnspan = 2)

    xscrollbar = Scrollbar(Tops, orient=HORIZONTAL)
    xscrollbar.grid(row=25, column=0, sticky=N+S+E+W,columnspan = 5)
    
    yscrollbar = Scrollbar(Tops)
    yscrollbar.grid(row=5, column=5, sticky=N+S+E+W,columnspan = 5)
    
    outputText = Text(Tops, wrap=NONE,height = 12, width = 90,font = ("Arial",16),bd = 5, 
                      xscrollcommand=xscrollbar.set,yscrollcommand=yscrollbar.set)      
    outputText.grid(row = 5,column = 0,columnspan = 5)
    outputText.configure(state ='disabled')
       
    root.mainloop()
