# Sequencing Objects 

In this module, we looked at sequences to determine their length, how many of each nucelotide, and even added and subtracted sequences. We were also able to translte vertain strands of sequences and use reverse complement. 
from Bio.Seq import Seq
```


```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq): 
    print("%i %s" % (index,letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# we can also print the length of each sequence 
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
Seq("AAAA").count("AA")

# we can also count how many times a certain sequence occurs
```




    2




```python
my_seq = Seq("GGATTCGATCACATGCCAGGAT")
```


```python
# we can count how long a sequence is
len(my_seq)
```




    22




```python
# we can also count how many of each letter in a seqeunce 
my_seq.count("G")
```




    6




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)

#we can determine what percentage of our sequence is what nucleotides
```




    50.0




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GGATTCGATCACATGCCAGGAT")
```


```python
gc_fraction(my_seq)
```




    0.5




```python
my_seq[4:12]
```




    Seq('TCGATCAC')




```python
my_seq[4:12]
```




    Seq('TCGATCAC')




```python
my_seq[0::3]
```




    Seq('GTGCACGT')




```python
my_seq[1::3]
```




    Seq('GTAATCG')




```python
my_seq[2:3]
```




    Seq('A')




```python
my_seq[::-1]
```




    Seq('TAGGACCGTACACTAGCTTAGG')




```python
str(my_seq)
```




    'GGATTCGATCACATGCCAGGAT'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GGATTCGATCACATGCCAGGAT
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")

```


```python
seq1 + seq2

# we can also add sequences
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")


```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCTATATAGGATCGAAAATCGC")
```


```python
my_seq.complement()
```




    Seq('CTAGCTACCCGATATATCCTAGCTTTTAGCG')




```python
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGCCCATCGATC')




```python
protein_seq = Seq("EVNARK")
protein_seq.complement()
```




    Seq('EBNTYM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna = coding_dna.reverse_complement().transcribe()
```


```python
template_dna
```




    Seq('CUAUCGGGCACCCUUUCAGCGGCCCAUUACAAUGGCCAU')



messenger_rna.back_transcribe()


```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
coding_dna.translate(table ="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table= 2)
```




    Seq('MAIVMGRWKGAR*')




```python
#we can  stop the sequence at the stop codon wihtout including it

coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table =2, to_stop = True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")

#we can also customize the stop symbol
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIGTRTFPGSGRSHGSTGCGNYVSPVNRWHLHGPPPPPRHHKKAPHDHH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIGTRTFPGSGRSHGSTGCGNYVSPVNRWHLHGPPPPPRHHKKAPHDHH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)

#we can change valine to methianine
```




    Seq('MKKMQSIGTRTFPGSGRSHGSTGCGNYVSPVNRWHLHGPPPPPRHHKKAPHDHH...HHR')



from Bio.Data import CodonTable 


```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGC"
```


```python
reverse_complement(my_string)
```




    'GCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGC'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGC'




```python
translate(my_string)
```

    /home/student/anaconda3/lib/python3.7/site-packages/Bio/Seq.py:2808: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
      BiopythonWarning,





    'AVMGRWKGGR'




```python

```


```python

```


```python

```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq): 
    print("%i %s" % (index,letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# we can also print the length of each sequence 
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
Seq("AAAA").count("AA")

# we can also count how many times a certain sequence occurs
```




    2




```python
my_seq = Seq("GGATTCGATCACATGCCAGGAT")
```


```python
# we can count how long a sequence is
len(my_seq)
```




    22




```python
# we can also count how many of each letter in a seqeunce 
my_seq.count("G")
```




    6




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)

#we can determine what percentage of our sequence is what nucleotides
```




    50.0




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GGATTCGATCACATGCCAGGAT")
```


```python
gc_fraction(my_seq)
```




    0.5




```python
my_seq[4:12]
```




    Seq('TCGATCAC')




```python
my_seq[4:12]
```




    Seq('TCGATCAC')




```python
my_seq[0::3]
```




    Seq('GTGCACGT')




```python
my_seq[1::3]
```




    Seq('GTAATCG')




```python
my_seq[2:3]
```




    Seq('A')




```python
my_seq[::-1]
```




    Seq('TAGGACCGTACACTAGCTTAGG')




```python
str(my_seq)
```




    'GGATTCGATCACATGCCAGGAT'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GGATTCGATCACATGCCAGGAT
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")

```


```python
seq1 + seq2

# we can also add sequences
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")


```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCTATATAGGATCGAAAATCGC")
```


```python
my_seq.complement()
```




    Seq('CTAGCTACCCGATATATCCTAGCTTTTAGCG')




```python
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGCCCATCGATC')




```python
protein_seq = Seq("EVNARK")
protein_seq.complement()
```




    Seq('EBNTYM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna = coding_dna.reverse_complement().transcribe()
```


```python
template_dna
```




    Seq('CUAUCGGGCACCCUUUCAGCGGCCCAUUACAAUGGCCAU')



messenger_rna.back_transcribe()


```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
coding_dna.translate(table ="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table= 2)
```




    Seq('MAIVMGRWKGAR*')




```python
#we can  stop the sequence at the stop codon wihtout including it

coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table =2, to_stop = True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")

#we can also customize the stop symbol
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIGTRTFPGSGRSHGSTGCGNYVSPVNRWHLHGPPPPPRHHKKAPHDHH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIGTRTFPGSGRSHGSTGCGNYVSPVNRWHLHGPPPPPRHHKKAPHDHH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)

#we can change valine to methianine
```




    Seq('MKKMQSIGTRTFPGSGRSHGSTGCGNYVSPVNRWHLHGPPPPPRHHKKAPHDHH...HHR')



from Bio.Data import CodonTable 


```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGC"
```


```python
reverse_complement(my_string)
```




    'GCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGC'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGC'




```python
translate(my_string)
```

    /home/student/anaconda3/lib/python3.7/site-packages/Bio/Seq.py:2808: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
      BiopythonWarning,





    'AVMGRWKGGR'




```python

```

# Sequence Annotations 

In this module, we looked at simple sequences and learned how to record them as well as examine their features and descriptions. 

from Bio.SeqRecord import SeqRecord
```


```python
help(SeqRecord)
```

    Help on class SeqRecord in module Bio.SeqRecord:
    
    class SeqRecord(builtins.object)
     |  SeqRecord(seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |  
     |  A SeqRecord object holds a sequence and information about it.
     |  
     |  Main attributes:
     |   - id          - Identifier such as a locus tag (string)
     |   - seq         - The sequence itself (Seq object or similar)
     |  
     |  Additional attributes:
     |   - name        - Sequence name, e.g. gene name (string)
     |   - description - Additional text (string)
     |   - dbxrefs     - List of database cross references (list of strings)
     |   - features    - Any (sub)features defined (list of SeqFeature objects)
     |   - annotations - Further information about the whole sequence (dictionary).
     |     Most entries are strings, or lists of strings.
     |   - letter_annotations - Per letter/symbol annotation (restricted
     |     dictionary). This holds Python sequences (lists, strings
     |     or tuples) whose length matches that of the sequence.
     |     A typical use would be to hold a list of integers
     |     representing sequencing quality scores, or a string
     |     representing the secondary structure.
     |  
     |  You will typically use Bio.SeqIO to read in sequences from files as
     |  SeqRecord objects.  However, you may want to create your own SeqRecord
     |  objects directly (see the __init__ method for further details):
     |  
     |  >>> from Bio.Seq import Seq
     |  >>> from Bio.SeqRecord import SeqRecord
     |  >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |  ...                    id="YP_025292.1", name="HokC",
     |  ...                    description="toxic membrane protein")
     |  >>> print(record)
     |  ID: YP_025292.1
     |  Name: HokC
     |  Description: toxic membrane protein
     |  Number of features: 0
     |  Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |  
     |  If you want to save SeqRecord objects to a sequence file, use Bio.SeqIO
     |  for this.  For the special case where you want the SeqRecord turned into
     |  a string in a particular file format there is a format method which uses
     |  Bio.SeqIO internally:
     |  
     |  >>> print(record.format("fasta"))
     |  >YP_025292.1 toxic membrane protein
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  <BLANKLINE>
     |  
     |  You can also do things like slicing a SeqRecord, checking its length, etc
     |  
     |  >>> len(record)
     |  44
     |  >>> edited = record[:10] + record[11:]
     |  >>> print(edited.seq)
     |  MKQHKAMIVAIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  >>> print(record.seq)
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  
     |  Methods defined here:
     |  
     |  __add__(self, other)
     |      Add another sequence or string to this sequence.
     |      
     |      The other sequence can be a SeqRecord object, a Seq object (or
     |      similar, e.g. a MutableSeq) or a plain Python string. If you add
     |      a plain string or a Seq (like) object, the new SeqRecord will simply
     |      have this appended to the existing data. However, any per letter
     |      annotation will be lost:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = record + "ACT"
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNACT
     |      >>> print(list(new.letter_annotations))
     |      []
     |      
     |      The new record will attempt to combine the annotation, but for any
     |      ambiguities (e.g. different names) it defaults to omitting that
     |      annotation.
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      
     |      Now let's cut the plasmid into two pieces, and join them back up the
     |      other way round (i.e. shift the starting point on this plasmid, have
     |      a look at the annotated features in the original file to see why this
     |      particular split point might make sense):
     |      
     |      >>> left = plasmid[:3765]
     |      >>> right = plasmid[3765:]
     |      >>> new = right + left
     |      >>> print("%s %i" % (new.id, len(new)))
     |      pBAD30 4923
     |      >>> str(new.seq) == str(right.seq + left.seq)
     |      True
     |      >>> len(new.features) == len(left.features) + len(right.features)
     |      True
     |      
     |      When we add the left and right SeqRecord objects, their annotation
     |      is all consistent, so it is all conserved in the new SeqRecord:
     |      
     |      >>> new.id == left.id == right.id == plasmid.id
     |      True
     |      >>> new.name == left.name == right.name == plasmid.name
     |      True
     |      >>> new.description == plasmid.description
     |      True
     |      >>> new.annotations == left.annotations == right.annotations
     |      True
     |      >>> new.letter_annotations == plasmid.letter_annotations
     |      True
     |      >>> new.dbxrefs == left.dbxrefs == right.dbxrefs
     |      True
     |      
     |      However, we should point out that when we sliced the SeqRecord,
     |      any annotations dictionary or dbxrefs list entries were lost.
     |      You can explicitly copy them like this:
     |      
     |      >>> new.annotations = plasmid.annotations.copy()
     |      >>> new.dbxrefs = plasmid.dbxrefs[:]
     |  
     |  __bool__(self)
     |      Boolean value of an instance of this class (True).
     |      
     |      This behaviour is for backwards compatibility, since until the
     |      __len__ method was added, a SeqRecord always evaluated as True.
     |      
     |      Note that in comparison, a Seq object will evaluate to False if it
     |      has a zero length sequence.
     |      
     |      WARNING: The SeqRecord may in future evaluate to False when its
     |      sequence is of zero length (in order to better match the Seq
     |      object behaviour)!
     |  
     |  __bytes__(self)
     |  
     |  __contains__(self, char)
     |      Implement the 'in' keyword, searches the sequence.
     |      
     |      e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> "GAATTC" in record
     |      False
     |      >>> "AAA" in record
     |      True
     |      
     |      This essentially acts as a proxy for using "in" on the sequence:
     |      
     |      >>> "GAATTC" in record.seq
     |      False
     |      >>> "AAA" in record.seq
     |      True
     |      
     |      Note that you can also use Seq objects as the query,
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> Seq("AAA") in record
     |      True
     |      
     |      See also the Seq object's __contains__ method.
     |  
     |  __eq__(self, other)
     |      Define the equal-to operand (not implemented).
     |  
     |  __format__(self, format_spec)
     |      Return the record as a string in the specified file format.
     |      
     |      This method supports the Python format() function and f-strings.
     |      The format_spec should be a lower case string supported by
     |      Bio.SeqIO as a text output file format. Requesting a binary file
     |      format raises a ValueError. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      ...
     |      >>> format(record, "fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(f"Here is {record.id} in FASTA format:\n{record:fasta}")
     |      Here is YP_025292.1 in FASTA format:
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      See also the SeqRecord's format() method.
     |  
     |  __ge__(self, other)
     |      Define the greater-than-or-equal-to operand (not implemented).
     |  
     |  __getitem__(self, index)
     |      Return a sub-sequence or an individual letter.
     |      
     |      Slicing, e.g. my_record[5:10], returns a new SeqRecord for
     |      that sub-sequence with some annotation preserved as follows:
     |      
     |      * The name, id and description are kept as-is.
     |      * Any per-letter-annotations are sliced to match the requested
     |        sub-sequence.
     |      * Unless a stride is used, all those features which fall fully
     |        within the subsequence are included (with their locations
     |        adjusted accordingly). If you want to preserve any truncated
     |        features (e.g. GenBank/EMBL source features), you must
     |        explicitly add them to the new SeqRecord yourself.
     |      * With the exception of any molecule type, the annotations
     |        dictionary and the dbxrefs list are not used for the new
     |        SeqRecord, as in general they may not apply to the
     |        subsequence. If you want to preserve them, you must explicitly
     |        copy them to the new SeqRecord yourself.
     |      
     |      Using an integer index, e.g. my_record[5] is shorthand for
     |      extracting that letter from the sequence, my_record.seq[5].
     |      
     |      For example, consider this short protein and its secondary
     |      structure as encoded by the PDB (e.g. H for alpha helices),
     |      plus a simple feature for its histidine self phosphorylation
     |      site:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
     |      >>> rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT"
     |      ...                     "EMMSEQDGYLAESINKDIEECNAIIEQFIDYLR"),
     |      ...                 id="1JOY", name="EnvZ",
     |      ...                 description="Homodimeric domain of EnvZ from E. coli")
     |      >>> rec.letter_annotations["secondary_structure"] = "  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  "
     |      >>> rec.features.append(SeqFeature(SimpleLocation(20, 21),
     |      ...                     type = "Site"))
     |      
     |      Now let's have a quick look at the full record,
     |      
     |      >>> print(rec)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      >>> rec.letter_annotations["secondary_structure"]
     |      '  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  '
     |      >>> print(rec.features[0].location)
     |      [20:21]
     |      
     |      Now let's take a sub sequence, here chosen as the first (fractured)
     |      alpha helix which includes the histidine phosphorylation site:
     |      
     |      >>> sub = rec[11:41]
     |      >>> print(sub)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('RTLLMAGVSHDLRTPLTRIRLATEMMSEQD')
     |      >>> sub.letter_annotations["secondary_structure"]
     |      'HHHHHTTTHHHHHHHHHHHHHHHHHHHHHH'
     |      >>> print(sub.features[0].location)
     |      [9:10]
     |      
     |      You can also of course omit the start or end values, for
     |      example to get the first ten letters only:
     |      
     |      >>> print(rec[:10])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLAD')
     |      
     |      Or for the last ten letters:
     |      
     |      >>> print(rec[-10:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('IIEQFIDYLR')
     |      
     |      If you omit both, then you get a copy of the original record (although
     |      lacking the annotations and dbxrefs):
     |      
     |      >>> print(rec[:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      
     |      Finally, indexing with a simple integer is shorthand for pulling out
     |      that letter from the sequence directly:
     |      
     |      >>> rec[5]
     |      'K'
     |      >>> rec.seq[5]
     |      'K'
     |  
     |  __gt__(self, other)
     |      Define the greater-than operand (not implemented).
     |  
     |  __init__(self, seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |      Create a SeqRecord.
     |      
     |      Arguments:
     |       - seq         - Sequence, required (Seq or MutableSeq)
     |       - id          - Sequence identifier, recommended (string)
     |       - name        - Sequence name, optional (string)
     |       - description - Sequence description, optional (string)
     |       - dbxrefs     - Database cross references, optional (list of strings)
     |       - features    - Any (sub)features, optional (list of SeqFeature objects)
     |       - annotations - Dictionary of annotations for the whole sequence
     |       - letter_annotations - Dictionary of per-letter-annotations, values
     |         should be strings, list or tuples of the same length as the full
     |         sequence.
     |      
     |      You will typically use Bio.SeqIO to read in sequences from files as
     |      SeqRecord objects.  However, you may want to create your own SeqRecord
     |      objects directly.
     |      
     |      Note that while an id is optional, we strongly recommend you supply a
     |      unique id string for each record.  This is especially important
     |      if you wish to write your sequences to a file.
     |      
     |      You can create a 'blank' SeqRecord object, and then populate the
     |      attributes later.
     |  
     |  __iter__(self)
     |      Iterate over the letters in the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a protein FASTA file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/loveliesbleeding.pro", "fasta")
     |      >>> for amino in record:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      This is just a shortcut for iterating over the sequence directly:
     |      
     |      >>> for amino in record.seq:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      Note that this does not facilitate iteration together with any
     |      per-letter-annotation.  However, you can achieve that using the
     |      python zip function on the record (or its sequence) and the relevant
     |      per-letter-annotation:
     |      
     |      >>> from Bio import SeqIO
     |      >>> rec = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(rec.letter_annotations))
     |      ['solexa_quality']
     |      >>> for nuc, qual in zip(rec, rec.letter_annotations["solexa_quality"]):
     |      ...     if qual > 35:
     |      ...         print("%s %i" % (nuc, qual))
     |      A 40
     |      C 39
     |      G 38
     |      T 37
     |      A 36
     |      
     |      You may agree that using zip(rec.seq, ...) is more explicit than using
     |      zip(rec, ...) as shown above.
     |  
     |  __le__(self, other)
     |      Define the less-than-or-equal-to operand (not implemented).
     |  
     |  __len__(self)
     |      Return the length of the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a FASTA nucleotide file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> len(record)
     |      309
     |      >>> len(record.seq)
     |      309
     |  
     |  __lt__(self, other)
     |      Define the less-than operand (not implemented).
     |  
     |  __ne__(self, other)
     |      Define the not-equal-to operand (not implemented).
     |  
     |  __radd__(self, other)
     |      Add another sequence or string to this sequence (from the left).
     |      
     |      This method handles adding a Seq object (or similar, e.g. MutableSeq)
     |      or a plain Python string (on the left) to a SeqRecord (on the right).
     |      See the __add__ method for more details, but for example:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = "ACT" + record
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(new.letter_annotations))
     |      []
     |  
     |  __repr__(self)
     |      Return a concise summary of the record for debugging (string).
     |      
     |      The python built in function repr works by calling the object's __repr__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT"
     |      ...                     "GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ"
     |      ...                     "SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ"
     |      ...                     "PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF"),
     |      ...                 id="NP_418483.1", name="b4059",
     |      ...                 description="ssDNA-binding protein",
     |      ...                 dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])
     |      >>> print(repr(rec))
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      At the python prompt you can also use this shorthand:
     |      
     |      >>> rec
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      Note that long sequences are shown truncated. Also note that any
     |      annotations, letter_annotations and features are not shown (as they
     |      would lead to a very long string).
     |  
     |  __str__(self)
     |      Return a human readable summary of the record and its annotation (string).
     |      
     |      The python built in function str works by calling the object's __str__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein, small")
     |      >>> print(str(record))
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      In this example you don't actually need to call str explicitly, as the
     |      print command does this automatically:
     |      
     |      >>> print(record)
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      Note that long sequences are shown truncated.
     |  
     |  count(self, sub, start=None, end=None)
     |      Return the number of non-overlapping occurrences of sub in seq[start:end].
     |      
     |      Optional arguments start and end are interpreted as in slice notation.
     |      This method behaves as the count method of Python strings.
     |  
     |  format(self, format)
     |      Return the record as a string in the specified file format.
     |      
     |      The format should be a lower case string supported as an output
     |      format by Bio.SeqIO, which is used to turn the SeqRecord into a
     |      string.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      >>> record.format("fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(record.format("fasta"))
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      The Python print function automatically appends a new line, meaning
     |      in this example a blank line is shown.  If you look at the string
     |      representation you can see there is a trailing new line (shown as
     |      slash n) which is important when writing to a file or if
     |      concatenating multiple sequence strings together.
     |      
     |      Note that this method will NOT work on every possible file format
     |      supported by Bio.SeqIO (e.g. some are for multiple sequences only,
     |      and binary formats are not supported).
     |  
     |  islower(self)
     |      Return True if all ASCII characters in the record's sequence are lowercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  isupper(self)
     |      Return True if all ASCII characters in the record's sequence are uppercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  lower(self)
     |      Return a copy of the record with a lower case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/aster.pro", "fasta")
     |      >>> print(record.format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVG
     |      VTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI
     |      <BLANKLINE>
     |      >>> print(record.lower().format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      gghvnpavtfgafvggnitllrgivyiiaqllgstvaclllkfvtndmavgvfslsagvg
     |      vtnalvfeivmtfglvytvyataidpkkgslgtiapiaigfivgani
     |      <BLANKLINE>
     |      
     |      To take a more annotation rich example,
     |      
     |      >>> from Bio import SeqIO
     |      >>> old = SeqIO.read("EMBL/TRBG361.embl", "embl")
     |      >>> len(old.features)
     |      3
     |      >>> new = old.lower()
     |      >>> len(old.features) == len(new.features)
     |      True
     |      >>> old.annotations["organism"] == new.annotations["organism"]
     |      True
     |      >>> old.dbxrefs == new.dbxrefs
     |      True
     |  
     |  reverse_complement(self, id=False, name=False, description=False, features=True, annotations=False, letter_annotations=True, dbxrefs=False)
     |      Return new SeqRecord with reverse complement sequence.
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the reversed sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or True to keep that of the parent, or False to
     |      omit them. The default is to keep the original features (with the
     |      strand and locations adjusted).
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent,
     |      or False to omit them. The default is to keep the original
     |      annotations (with the letter annotations reversed).
     |      
     |      To show what happens to the pre-letter annotations, consider an
     |      example Solexa variant FASTQ file with a single entry, which we'll
     |      read in as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Now take the reverse complement, here we explicitly give a new
     |      identifier (the old identifier with a suffix):
     |      
     |      >>> rc_record = record.reverse_complement(id=record.id + "_rc")
     |      >>> print("%s %s" % (rc_record.id, rc_record.seq))
     |      slxa_0001_1_0001_01_rc NNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
     |      
     |      Notice that the per-letter-annotations have also been reversed,
     |      although this may not be appropriate for all cases.
     |      
     |      >>> print(rc_record.letter_annotations["solexa_quality"])
     |      [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
     |      
     |      Now for the features, we need a different example. Parsing a GenBank
     |      file is probably the easiest way to get an nice example with features
     |      in it...
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      >>> plasmid.seq
     |      Seq('GCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGA...ATG')
     |      >>> len(plasmid.features)
     |      13
     |      
     |      Now, let's take the reverse complement of this whole plasmid:
     |      
     |      >>> rc_plasmid = plasmid.reverse_complement(id=plasmid.id+"_rc")
     |      >>> print("%s %i" % (rc_plasmid.id, len(rc_plasmid)))
     |      pBAD30_rc 4923
     |      >>> rc_plasmid.seq
     |      Seq('CATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCA...AGC')
     |      >>> len(rc_plasmid.features)
     |      13
     |      
     |      Let's compare the first CDS feature - it has gone from being the
     |      second feature (index 1) to the second last feature (index -2), its
     |      strand has changed, and the location switched round.
     |      
     |      >>> print(plasmid.features[1])
     |      type: CDS
     |      location: [1081:1960](-)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      >>> print(rc_plasmid.features[-2])
     |      type: CDS
     |      location: [2963:3842](+)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      
     |      You can check this new location, based on the length of the plasmid:
     |      
     |      >>> len(plasmid) - 1081
     |      3842
     |      >>> len(plasmid) - 1960
     |      2963
     |      
     |      Note that if the SeqFeature annotation includes any strand specific
     |      information (e.g. base changes for a SNP), this information is not
     |      amended, and would need correction after the reverse complement.
     |      
     |      Note trying to reverse complement a protein SeqRecord raises an
     |      exception:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> protein_rec = SeqRecord(Seq("MAIVMGR"), id="Test",
     |      ...                         annotations={"molecule_type": "protein"})
     |      >>> protein_rec.reverse_complement()
     |      Traceback (most recent call last):
     |         ...
     |      ValueError: Proteins do not have complements!
     |      
     |      If you have RNA without any U bases, it must be annotated as RNA
     |      otherwise it will be treated as DNA by default with A mapped to T:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rna1 = SeqRecord(Seq("ACG"), id="Test")
     |      >>> rna2 = SeqRecord(Seq("ACG"), id="Test", annotations={"molecule_type": "RNA"})
     |      >>> print(rna1.reverse_complement(id="RC", description="unk").format("fasta"))
     |      >RC unk
     |      CGT
     |      <BLANKLINE>
     |      >>> print(rna2.reverse_complement(id="RC", description="RNA").format("fasta"))
     |      >RC RNA
     |      CGU
     |      <BLANKLINE>
     |      
     |      Also note you can reverse complement a SeqRecord using a MutableSeq:
     |      
     |      >>> from Bio.Seq import MutableSeq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(MutableSeq("ACGT"), id="Test")
     |      >>> rec.seq[0] = "T"
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      Test TCGT
     |      >>> rc = rec.reverse_complement(id=True)
     |      >>> print("%s %s" % (rc.id, rc.seq))
     |      Test ACGA
     |  
     |  translate(self, table='Standard', stop_symbol='*', to_stop=False, cds=False, gap=None, id=False, name=False, description=False, features=False, annotations=False, letter_annotations=False, dbxrefs=False)
     |      Return new SeqRecord with translated sequence.
     |      
     |      This calls the record's .seq.translate() method (which describes
     |      the translation related arguments, like table for the genetic code),
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the translated sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or False (default) to omit them.
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent
     |      (annotations only), or False (default) to omit them.
     |      
     |      e.g. Loading a FASTA gene and translating it,
     |      
     |      >>> from Bio import SeqIO
     |      >>> gene_record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> print(gene_record.format("fasta"))
     |      >gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds
     |      CAGGCTGCGCGGTTTCTATTTATGAAGAACAAGGTCCGTATGATAGTTGATTGTCATGCA
     |      AAACATGTGAAGGTTCTTCAAGACGAAAAACTCCCATTTGATTTGACTCTGTGCGGTTCG
     |      ACCTTAAGAGCTCCACATAGTTGCCATTTGCAGTACATGGCTAACATGGATTCAATTGCT
     |      TCATTGGTTATGGCAGTGGTCGTCAATGACAGCGATGAAGATGGAGATAGCCGTGACGCA
     |      GTTCTACCACAAAAGAAAAAGAGACTTTGGGGTTTGGTAGTTTGTCATAACACTACTCCG
     |      AGGTTTGTT
     |      <BLANKLINE>
     |      
     |      And now translating the record, specifying the new ID and description:
     |      
     |      >>> protein_record = gene_record.translate(table=11,
     |      ...                                        id="phya",
     |      ...                                        description="translation")
     |      >>> print(protein_record.format("fasta"))
     |      >phya translation
     |      QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
     |      SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV
     |      <BLANKLINE>
     |  
     |  upper(self)
     |      Return a copy of the record with an upper case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("acgtACGT"), id="Test",
     |      ...                    description = "Made up for this example")
     |      >>> record.letter_annotations["phred_quality"] = [1, 2, 3, 4, 5, 6, 7, 8]
     |      >>> print(record.upper().format("fastq"))
     |      @Test Made up for this example
     |      ACGTACGT
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |      
     |      Naturally, there is a matching lower method:
     |      
     |      >>> print(record.lower().format("fastq"))
     |      @Test Made up for this example
     |      acgtacgt
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  letter_annotations
     |      Dictionary of per-letter-annotation for the sequence.
     |      
     |      For example, this can hold quality scores used in FASTQ or QUAL files.
     |      Consider this example using Bio.SeqIO to read in an example Solexa
     |      variant FASTQ file as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      The letter_annotations get sliced automatically if you slice the
     |      parent SeqRecord, for example taking the last ten bases:
     |      
     |      >>> sub_record = record[-10:]
     |      >>> print("%s %s" % (sub_record.id, sub_record.seq))
     |      slxa_0001_1_0001_01 ACGTNNNNNN
     |      >>> print(sub_record.letter_annotations["solexa_quality"])
     |      [4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Any python sequence (i.e. list, tuple or string) can be recorded in
     |      the SeqRecord's letter_annotations dictionary as long as the length
     |      matches that of the SeqRecord's sequence.  e.g.
     |      
     |      >>> len(sub_record.letter_annotations)
     |      1
     |      >>> sub_record.letter_annotations["dummy"] = "abcdefghij"
     |      >>> len(sub_record.letter_annotations)
     |      2
     |      
     |      You can delete entries from the letter_annotations dictionary as usual:
     |      
     |      >>> del sub_record.letter_annotations["solexa_quality"]
     |      >>> sub_record.letter_annotations
     |      {'dummy': 'abcdefghij'}
     |      
     |      You can completely clear the dictionary easily as follows:
     |      
     |      >>> sub_record.letter_annotations = {}
     |      >>> sub_record.letter_annotations
     |      {}
     |      
     |      Note that if replacing the record's sequence with a sequence of a
     |      different length you must first clear the letter_annotations dict.
     |  
     |  seq
     |      The sequence itself, as a Seq or MutableSeq object.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __hash__ = None
    



```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB computational biology class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB computational biology class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB computational biology class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"]= "None. This is just an example."
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example.



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB computational biology class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_quality"]=[40,40,38,30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
#https://raw.githubusercontent.com/biopyhton/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
#https://raw.githubusercontent.com/biopyhton/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left =8, right =9)
```


```python

```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
...                    "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
...                    "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
...                    "SSAC") , 
                              id="gi|14150838|gb|AAK54648.1|AF376133_1",
...                description= "chalcone synthase [Cucumis sativus]", 
...  )
```


```python
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
len(sub_record.features)
```




    2




```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description = 'Yersinia pestis biovr Microtus str. 91001 plasmid pPCP1, complete sequence'
```


```python
sub_record.description
```




    'Yersinia pestis biovr Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
print(sub_record.format("genbank"[:200] + "..."))
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-184-1f20f57b3390> in <module>
    ----> 1 print(sub_record.format("genbank"[:200] + "..."))
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqRecord.py in format(self, format)
        726         """
        727         # See also the __format__ method
    --> 728         return self.__format__(format)
        729 
        730     def __format__(self, format_spec):


    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqRecord.py in __format__(self, format_spec)
        764         handle = StringIO()
        765         try:
    --> 766             SeqIO.write(self, handle, format_spec)
        767         except StreamModeError:
        768             raise ValueError(


    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/__init__.py in write(sequences, handle, format)
        546         raise ValueError(f"Reading format '{format}' is supported, but not writing")
        547 
    --> 548     raise ValueError(f"Unknown format '{format}'")
        549 
        550 


    ValueError: Unknown format 'genbank...'



```python
record =SeqIO.read("NC_005816.gb.txt" , "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i " % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13 



```python
rc = record.reverse_complement(id= "Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



```python

```
# Sequence Input and Output 

In this module, we looked at fasta files and printed sequences. We also incorporated annotations, descriptions, and records. 

#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk
```


```python
from Bio import SeqIO
```


```python
help(SeqIO)
```

    Help on package Bio.SeqIO in Bio:
    
    NAME
        Bio.SeqIO - Sequence input/output as SeqRecord objects.
    
    DESCRIPTION
        Bio.SeqIO is also documented at SeqIO_ and by a whole chapter in our tutorial:
        
          - `HTML Tutorial`_
          - `PDF Tutorial`_
        
        .. _SeqIO: http://biopython.org/wiki/SeqIO
        .. _`HTML Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.html
        .. _`PDF Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.pdf
        
        Input
        -----
        The main function is Bio.SeqIO.parse(...) which takes an input file handle
        (or in recent versions of Biopython alternatively a filename as a string),
        and format string.  This returns an iterator giving SeqRecord objects:
        
        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("Fasta/f002", "fasta"):
        ...     print("%s %i" % (record.id, len(record)))
        gi|1348912|gb|G26680|G26680 633
        gi|1348917|gb|G26685|G26685 413
        gi|1592936|gb|G29385|G29385 471
        
        Note that the parse() function will invoke the relevant parser for the
        format with its default settings.  You may want more control, in which case
        you need to create a format specific sequence iterator directly.
        
        Some of these parsers are wrappers around low-level parsers which build up
        SeqRecord objects for the consistent SeqIO interface. In cases where the
        run-time is critical, such as large FASTA or FASTQ files, calling these
        underlying parsers will be much faster - in this case these generator
        functions which return tuples of strings:
        
        >>> from Bio.SeqIO.FastaIO import SimpleFastaParser
        >>> from Bio.SeqIO.QualityIO import FastqGeneralIterator
        
        
        Input - Single Records
        ----------------------
        If you expect your file to contain one-and-only-one record, then we provide
        the following 'helper' function which will return a single SeqRecord, or
        raise an exception if there are no records or more than one record:
        
        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Fasta/f001", "fasta")
        >>> print("%s %i" % (record.id, len(record)))
        gi|3318709|pdb|1A91| 79
        
        This style is useful when you expect a single record only (and would
        consider multiple records an error).  For example, when dealing with GenBank
        files for bacterial genomes or chromosomes, there is normally only a single
        record.  Alternatively, use this with a handle when downloading a single
        record from the internet.
        
        However, if you just want the first record from a file containing multiple
        record, use the next() function on the iterator:
        
        >>> from Bio import SeqIO
        >>> record = next(SeqIO.parse("Fasta/f002", "fasta"))
        >>> print("%s %i" % (record.id, len(record)))
        gi|1348912|gb|G26680|G26680 633
        
        The above code will work as long as the file contains at least one record.
        Note that if there is more than one record, the remaining records will be
        silently ignored.
        
        
        Input - Multiple Records
        ------------------------
        For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records
        using a sequence iterator can save you a lot of memory (RAM).  There is
        less benefit for interlaced file formats (e.g. most multiple alignment file
        formats).  However, an iterator only lets you access the records one by one.
        
        If you want random access to the records by number, turn this into a list:
        
        >>> from Bio import SeqIO
        >>> records = list(SeqIO.parse("Fasta/f002", "fasta"))
        >>> len(records)
        3
        >>> print(records[1].id)
        gi|1348917|gb|G26685|G26685
        
        If you want random access to the records by a key such as the record id,
        turn the iterator into a dictionary:
        
        >>> from Bio import SeqIO
        >>> record_dict = SeqIO.to_dict(SeqIO.parse("Fasta/f002", "fasta"))
        >>> len(record_dict)
        3
        >>> print(len(record_dict["gi|1348917|gb|G26685|G26685"]))
        413
        
        However, using list() or the to_dict() function will load all the records
        into memory at once, and therefore is not possible on very large files.
        Instead, for *some* file formats Bio.SeqIO provides an indexing approach
        providing dictionary like access to any record. For example,
        
        >>> from Bio import SeqIO
        >>> record_dict = SeqIO.index("Fasta/f002", "fasta")
        >>> len(record_dict)
        3
        >>> print(len(record_dict["gi|1348917|gb|G26685|G26685"]))
        413
        >>> record_dict.close()
        
        Many but not all of the supported input file formats can be indexed like
        this. For example "fasta", "fastq", "qual" and even the binary format "sff"
        work, but alignment formats like "phylip", "clustalw" and "nexus" will not.
        
        In most cases you can also use SeqIO.index to get the record from the file
        as a raw string (not a SeqRecord). This can be useful for example to extract
        a sub-set of records from a file where SeqIO cannot output the file format
        (e.g. the plain text SwissProt format, "swiss") or where it is important to
        keep the output 100% identical to the input). For example,
        
        >>> from Bio import SeqIO
        >>> record_dict = SeqIO.index("Fasta/f002", "fasta")
        >>> len(record_dict)
        3
        >>> print(record_dict.get_raw("gi|1348917|gb|G26685|G26685").decode())
        >gi|1348917|gb|G26685|G26685 human STS STS_D11734.
        CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATCTTACTCTTTC
        ATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTTTACAGATGTGAAACTTTCAA
        GAGATTTACTGACTTTCCTAGAATAGTTTCTCTACTGGAAACCTGATGCTTTTATAAGCCATTGTGATTA
        GGATGACTGTTACAGGCTTAGCTTTGTGTGAAANCCAGTCACCTTTCTCCTAGGTAATGAGTAGTGCTGT
        TCATATTACTNTAAGTTCTATAGCATACTTGCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGG
        AATTGNTTTGCCCTTTTGGGTTTNTTNTTGGTAAANNNTTCCCGGGTGGGGGNGGTNNNGAAA
        <BLANKLINE>
        >>> print(record_dict["gi|1348917|gb|G26685|G26685"].format("fasta"))
        >gi|1348917|gb|G26685|G26685 human STS STS_D11734.
        CGGAGCCAGCGAGCATATGCTGCATGAGGACCTTTCTATCTTACATTATGGCTGGGAATC
        TTACTCTTTCATCTGATACCTTGTTCAGATTTCAAAATAGTTGTAGCCTTATCCTGGTTT
        TACAGATGTGAAACTTTCAAGAGATTTACTGACTTTCCTAGAATAGTTTCTCTACTGGAA
        ACCTGATGCTTTTATAAGCCATTGTGATTAGGATGACTGTTACAGGCTTAGCTTTGTGTG
        AAANCCAGTCACCTTTCTCCTAGGTAATGAGTAGTGCTGTTCATATTACTNTAAGTTCTA
        TAGCATACTTGCNATCCTTTANCCATGCTTATCATANGTACCATTTGAGGAATTGNTTTG
        CCCTTTTGGGTTTNTTNTTGGTAAANNNTTCCCGGGTGGGGGNGGTNNNGAAA
        <BLANKLINE>
        >>> record_dict.close()
        
        Here the original file and what Biopython would output differ in the line
        wrapping. Also note that the get_raw method will return a bytes object,
        hence the use of decode to turn it into a string.
        
        Also note that the get_raw method will preserve the newline endings. This
        example FASTQ file uses Unix style endings (b"\n" only),
        
        >>> from Bio import SeqIO
        >>> fastq_dict = SeqIO.index("Quality/example.fastq", "fastq")
        >>> len(fastq_dict)
        3
        >>> raw = fastq_dict.get_raw("EAS54_6_R1_2_1_540_792")
        >>> raw.count(b"\n")
        4
        >>> raw.count(b"\r\n")
        0
        >>> b"\r" in raw
        False
        >>> len(raw)
        78
        >>> fastq_dict.close()
        
        Here is the same file but using DOS/Windows new lines (b"\r\n" instead),
        
        >>> from Bio import SeqIO
        >>> fastq_dict = SeqIO.index("Quality/example_dos.fastq", "fastq")
        >>> len(fastq_dict)
        3
        >>> raw = fastq_dict.get_raw("EAS54_6_R1_2_1_540_792")
        >>> raw.count(b"\n")
        4
        >>> raw.count(b"\r\n")
        4
        >>> b"\r\n" in raw
        True
        >>> len(raw)
        82
        >>> fastq_dict.close()
        
        Because this uses two bytes for each new line, the file is longer than
        the Unix equivalent with only one byte.
        
        
        Input - Alignments
        ------------------
        You can read in alignment files as alignment objects using Bio.AlignIO.
        Alternatively, reading in an alignment file format via Bio.SeqIO will give
        you a SeqRecord for each row of each alignment:
        
        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("Clustalw/hedgehog.aln", "clustal"):
        ...     print("%s %i" % (record.id, len(record)))
        gi|167877390|gb|EDS40773.1| 447
        gi|167234445|ref|NP_001107837. 447
        gi|74100009|gb|AAZ99217.1| 447
        gi|13990994|dbj|BAA33523.2| 447
        gi|56122354|gb|AAV74328.1| 447
        
        
        Output
        ------
        Use the function Bio.SeqIO.write(...), which takes a complete set of
        SeqRecord objects (either as a list, or an iterator), an output file handle
        (or in recent versions of Biopython an output filename as a string) and of
        course the file format::
        
          from Bio import SeqIO
          records = ...
          SeqIO.write(records, "example.faa", "fasta")
        
        Or, using a handle::
        
            from Bio import SeqIO
            records = ...
            with open("example.faa", "w") as handle:
              SeqIO.write(records, handle, "fasta")
        
        You are expected to call this function once (with all your records) and if
        using a handle, make sure you close it to flush the data to the hard disk.
        
        
        Output - Advanced
        -----------------
        The effect of calling write() multiple times on a single file will vary
        depending on the file format, and is best avoided unless you have a strong
        reason to do so.
        
        If you give a filename, then each time you call write() the existing file
        will be overwritten. For sequential files formats (e.g. fasta, genbank) each
        "record block" holds a single sequence.  For these files it would probably
        be safe to call write() multiple times by re-using the same handle.
        
        However, trying this for certain alignment formats (e.g. phylip, clustal,
        stockholm) would have the effect of concatenating several multiple sequence
        alignments together.  Such files are created by the PHYLIP suite of programs
        for bootstrap analysis, but it is clearer to do this via Bio.AlignIO instead.
        
        Worse, many fileformats have an explicit header and/or footer structure
        (e.g. any XMl format, and most binary file formats like SFF). Here making
        multiple calls to write() will result in an invalid file.
        
        
        Conversion
        ----------
        The Bio.SeqIO.convert(...) function allows an easy interface for simple
        file format conversions. Additionally, it may use file format specific
        optimisations so this should be the fastest way too.
        
        In general however, you can combine the Bio.SeqIO.parse(...) function with
        the Bio.SeqIO.write(...) function for sequence file conversion. Using
        generator expressions or generator functions provides a memory efficient way
        to perform filtering or other extra operations as part of the process.
        
        
        File Formats
        ------------
        When specifying the file format, use lowercase strings.  The same format
        names are also used in Bio.AlignIO and include the following:
        
            - abi     - Applied Biosystem's sequencing trace format
            - abi-trim - Same as "abi" but with quality trimming with Mott's algorithm
            - ace     - Reads the contig sequences from an ACE assembly file.
            - cif-atom - Uses Bio.PDB.MMCIFParser to determine the (partial) protein
              sequence as it appears in the structure based on the atomic coordinates.
            - cif-seqres - Reads a macromolecular Crystallographic Information File
              (mmCIF) file to determine the complete protein sequence as defined by the
              _pdbx_poly_seq_scheme records.
            - embl    - The EMBL flat file format. Uses Bio.GenBank internally.
            - fasta   - The generic sequence file format where each record starts with
              an identifier line starting with a ">" character, followed by
              lines of sequence.
            - fasta-2line - Stricter interpretation of the FASTA format using exactly
              two lines per record (no line wrapping).
            - fastq   - A "FASTA like" format used by Sanger which also stores PHRED
              sequence quality values (with an ASCII offset of 33).
            - fastq-sanger - An alias for "fastq" for consistency with BioPerl and EMBOSS
            - fastq-solexa - Original Solexa/Illumnia variant of the FASTQ format which
              encodes Solexa quality scores (not PHRED quality scores) with an
              ASCII offset of 64.
            - fastq-illumina - Solexa/Illumina 1.3 to 1.7 variant of the FASTQ format
              which encodes PHRED quality scores with an ASCII offset of 64
              (not 33). Note as of version 1.8 of the CASAVA pipeline Illumina
              will produce FASTQ files using the standard Sanger encoding.
            - gck     - Gene Construction Kit's format.
            - genbank - The GenBank or GenPept flat file format.
            - gb      - An alias for "genbank", for consistency with NCBI Entrez Utilities
            - ig      - The IntelliGenetics file format, apparently the same as the
              MASE alignment format.
            - imgt    - An EMBL like format from IMGT where the feature tables are more
              indented to allow for longer feature types.
            - nib     - UCSC's nib file format for nucleotide sequences, which uses one
              nibble (4 bits) to represent each nucleotide, and stores two nucleotides in
              one byte.
            - pdb-seqres -  Reads a Protein Data Bank (PDB) file to determine the
              complete protein sequence as it appears in the header (no dependencies).
            - pdb-atom - Uses Bio.PDB to determine the (partial) protein sequence as
              it appears in the structure based on the atom coordinate section of the
              file (requires NumPy for Bio.PDB).
            - phd     - Output from PHRED, used by PHRAP and CONSED for input.
            - pir     - A "FASTA like" format introduced by the National Biomedical
              Research Foundation (NBRF) for the Protein Information Resource
              (PIR) database, now part of UniProt.
            - seqxml  - SeqXML, simple XML format described in Schmitt et al (2011).
            - sff     - Standard Flowgram Format (SFF), typical output from Roche 454.
            - sff-trim - Standard Flowgram Format (SFF) with given trimming applied.
            - snapgene - SnapGene's native format.
            - swiss   - Plain text Swiss-Prot aka UniProt format.
            - tab     - Simple two column tab separated sequence files, where each
              line holds a record's identifier and sequence. For example,
              this is used as by Aligent's eArray software when saving
              microarray probes in a minimal tab delimited text file.
            - qual    - A "FASTA like" format holding PHRED quality values from
              sequencing DNA, but no actual sequences (usually provided
              in separate FASTA files).
            - uniprot-xml - The UniProt XML format (replacement for the SwissProt plain
              text format which we call "swiss")
            - xdna        - DNA Strider's and SerialCloner's native format.
        
        Note that while Bio.SeqIO can read all the above file formats, it cannot
        write to all of them.
        
        You can also use any file format supported by Bio.AlignIO, such as "nexus",
        "phylip" and "stockholm", which gives you access to the individual sequences
        making up each alignment as SeqRecords.
    
    PACKAGE CONTENTS
        AbiIO
        AceIO
        FastaIO
        GckIO
        IgIO
        InsdcIO
        Interfaces
        NibIO
        PdbIO
        PhdIO
        PirIO
        QualityIO
        SeqXmlIO
        SffIO
        SnapGeneIO
        SwissIO
        TabIO
        TwoBitIO
        UniprotIO
        XdnaIO
        _index
        _twoBitIO
    
    FUNCTIONS
        convert(in_file, in_format, out_file, out_format, molecule_type=None)
            Convert between two sequence file formats, return number of records.
            
            Arguments:
             - in_file - an input handle or filename
             - in_format - input file format, lower case string
             - out_file - an output handle or filename
             - out_format - output file format, lower case string
             - molecule_type - optional molecule type to apply, string containing
               "DNA", "RNA" or "protein".
            
            **NOTE** - If you provide an output filename, it will be opened which will
            overwrite any existing file without warning.
            
            The idea here is that while doing this will work::
            
                from Bio import SeqIO
                records = SeqIO.parse(in_handle, in_format)
                count = SeqIO.write(records, out_handle, out_format)
            
            it is shorter to write::
            
                from Bio import SeqIO
                count = SeqIO.convert(in_handle, in_format, out_handle, out_format)
            
            Also, Bio.SeqIO.convert is faster for some conversions as it can make some
            optimisations.
            
            For example, going from a filename to a handle:
            
            >>> from Bio import SeqIO
            >>> from io import StringIO
            >>> handle = StringIO("")
            >>> SeqIO.convert("Quality/example.fastq", "fastq", handle, "fasta")
            3
            >>> print(handle.getvalue())
            >EAS54_6_R1_2_1_413_324
            CCCTTCTTGTCTTCAGCGTTTCTCC
            >EAS54_6_R1_2_1_540_792
            TTGGCAGGCCAAGGCCGATGGATCA
            >EAS54_6_R1_2_1_443_348
            GTTGCTTCTGGCGTGGGTGGGGGGG
            <BLANKLINE>
            
            Note some formats like SeqXML require you to specify the molecule type
            when it cannot be determined by the parser:
            
            >>> from Bio import SeqIO
            >>> from io import BytesIO
            >>> handle = BytesIO()
            >>> SeqIO.convert("Quality/example.fastq", "fastq", handle, "seqxml", "DNA")
            3
        
        index(filename, format, alphabet=None, key_function=None)
            Indexes a sequence file and returns a dictionary like object.
            
            Arguments:
             - filename - string giving name of file to be indexed
             - format   - lower case string describing the file format
             - alphabet - no longer used, leave as None
             - key_function - Optional callback function which when given a
               SeqRecord identifier string should return a unique key for the
               dictionary.
            
            This indexing function will return a dictionary like object, giving the
            SeqRecord objects as values.
            
            As of Biopython 1.69, this will preserve the ordering of the records in
            file when iterating over the entries.
            
            >>> from Bio import SeqIO
            >>> records = SeqIO.index("Quality/example.fastq", "fastq")
            >>> len(records)
            3
            >>> list(records)  # make a list of the keys
            ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_540_792', 'EAS54_6_R1_2_1_443_348']
            >>> print(records["EAS54_6_R1_2_1_540_792"].format("fasta"))
            >EAS54_6_R1_2_1_540_792
            TTGGCAGGCCAAGGCCGATGGATCA
            <BLANKLINE>
            >>> "EAS54_6_R1_2_1_540_792" in records
            True
            >>> print(records.get("Missing", None))
            None
            >>> records.close()
            
            If the file is BGZF compressed, this is detected automatically. Ordinary
            GZIP files are not supported:
            
            >>> from Bio import SeqIO
            >>> records = SeqIO.index("Quality/example.fastq.bgz", "fastq")
            >>> len(records)
            3
            >>> print(records["EAS54_6_R1_2_1_540_792"].seq)
            TTGGCAGGCCAAGGCCGATGGATCA
            >>> records.close()
            
            When you call the index function, it will scan through the file, noting
            the location of each record. When you access a particular record via the
            dictionary methods, the code will jump to the appropriate part of the
            file and then parse that section into a SeqRecord.
            
            Note that not all the input formats supported by Bio.SeqIO can be used
            with this index function. It is designed to work only with sequential
            file formats (e.g. "fasta", "gb", "fastq") and is not suitable for any
            interlaced file format (e.g. alignment formats such as "clustal").
            
            For small files, it may be more efficient to use an in memory Python
            dictionary, e.g.
            
            >>> from Bio import SeqIO
            >>> records = SeqIO.to_dict(SeqIO.parse("Quality/example.fastq", "fastq"))
            >>> len(records)
            3
            >>> list(records)  # make a list of the keys
            ['EAS54_6_R1_2_1_413_324', 'EAS54_6_R1_2_1_540_792', 'EAS54_6_R1_2_1_443_348']
            >>> print(records["EAS54_6_R1_2_1_540_792"].format("fasta"))
            >EAS54_6_R1_2_1_540_792
            TTGGCAGGCCAAGGCCGATGGATCA
            <BLANKLINE>
            
            As with the to_dict() function, by default the id string of each record
            is used as the key. You can specify a callback function to transform
            this (the record identifier string) into your preferred key. For example:
            
            >>> from Bio import SeqIO
            >>> def make_tuple(identifier):
            ...     parts = identifier.split("_")
            ...     return int(parts[-2]), int(parts[-1])
            >>> records = SeqIO.index("Quality/example.fastq", "fastq",
            ...                       key_function=make_tuple)
            >>> len(records)
            3
            >>> list(records)  # make a list of the keys
            [(413, 324), (540, 792), (443, 348)]
            >>> print(records[(540, 792)].format("fasta"))
            >EAS54_6_R1_2_1_540_792
            TTGGCAGGCCAAGGCCGATGGATCA
            <BLANKLINE>
            >>> (540, 792) in records
            True
            >>> "EAS54_6_R1_2_1_540_792" in records
            False
            >>> print(records.get("Missing", None))
            None
            >>> records.close()
            
            Another common use case would be indexing an NCBI style FASTA file,
            where you might want to extract the GI number from the FASTA identifier
            to use as the dictionary key.
            
            Notice that unlike the to_dict() function, here the key_function does
            not get given the full SeqRecord to use to generate the key. Doing so
            would impose a severe performance penalty as it would require the file
            to be completely parsed while building the index. Right now this is
            usually avoided.
            
            See Also: Bio.SeqIO.index_db() and Bio.SeqIO.to_dict()
        
        index_db(index_filename, filenames=None, format=None, alphabet=None, key_function=None)
            Index several sequence files and return a dictionary like object.
            
            The index is stored in an SQLite database rather than in memory (as in the
            Bio.SeqIO.index(...) function).
            
            Arguments:
             - index_filename - Where to store the SQLite index
             - filenames - list of strings specifying file(s) to be indexed, or when
               indexing a single file this can be given as a string.
               (optional if reloading an existing index, but must match)
             - format   - lower case string describing the file format
               (optional if reloading an existing index, but must match)
             - alphabet - no longer used, leave as None.
             - key_function - Optional callback function which when given a
               SeqRecord identifier string should return a unique
               key for the dictionary.
            
            This indexing function will return a dictionary like object, giving the
            SeqRecord objects as values:
            
            >>> from Bio import SeqIO
            >>> files = ["GenBank/NC_000932.faa", "GenBank/NC_005816.faa"]
            >>> def get_gi(name):
            ...     parts = name.split("|")
            ...     i = parts.index("gi")
            ...     assert i != -1
            ...     return parts[i+1]
            >>> idx_name = ":memory:" #use an in memory SQLite DB for this test
            >>> records = SeqIO.index_db(idx_name, files, "fasta", key_function=get_gi)
            >>> len(records)
            95
            >>> records["7525076"].description
            'gi|7525076|ref|NP_051101.1| Ycf2 [Arabidopsis thaliana]'
            >>> records["45478717"].description
            'gi|45478717|ref|NP_995572.1| pesticin [Yersinia pestis biovar Microtus str. 91001]'
            >>> records.close()
            
            In this example the two files contain 85 and 10 records respectively.
            
            BGZF compressed files are supported, and detected automatically. Ordinary
            GZIP compressed files are not supported.
            
            See Also: Bio.SeqIO.index() and Bio.SeqIO.to_dict(), and the Python module
            glob which is useful for building lists of files.
        
        parse(handle, format, alphabet=None)
            Turn a sequence file into an iterator returning SeqRecords.
            
            Arguments:
             - handle   - handle to the file, or the filename as a string
               (note older versions of Biopython only took a handle).
             - format   - lower case string describing the file format.
             - alphabet - no longer used, should be None.
            
            Typical usage, opening a file to read in, and looping over the record(s):
            
            >>> from Bio import SeqIO
            >>> filename = "Fasta/sweetpea.nu"
            >>> for record in SeqIO.parse(filename, "fasta"):
            ...    print("ID %s" % record.id)
            ...    print("Sequence length %i" % len(record))
            ID gi|3176602|gb|U78617.1|LOU78617
            Sequence length 309
            
            For lazy-loading file formats such as twobit, for which the file contents
            is read on demand only, ensure that the file remains open while extracting
            sequence data.
            
            If you have a string 'data' containing the file contents, you must
            first turn this into a handle in order to parse it:
            
            >>> data = ">Alpha\nACCGGATGTA\n>Beta\nAGGCTCGGTTA\n"
            >>> from Bio import SeqIO
            >>> from io import StringIO
            >>> for record in SeqIO.parse(StringIO(data), "fasta"):
            ...     print("%s %s" % (record.id, record.seq))
            Alpha ACCGGATGTA
            Beta AGGCTCGGTTA
            
            Use the Bio.SeqIO.read(...) function when you expect a single record
            only.
        
        read(handle, format, alphabet=None)
            Turn a sequence file into a single SeqRecord.
            
            Arguments:
             - handle   - handle to the file, or the filename as a string
               (note older versions of Biopython only took a handle).
             - format   - string describing the file format.
             - alphabet - no longer used, should be None.
            
            This function is for use parsing sequence files containing
            exactly one record.  For example, reading a GenBank file:
            
            >>> from Bio import SeqIO
            >>> record = SeqIO.read("GenBank/arab1.gb", "genbank")
            >>> print("ID %s" % record.id)
            ID AC007323.5
            >>> print("Sequence length %i" % len(record))
            Sequence length 86436
            
            If the handle contains no records, or more than one record,
            an exception is raised.  For example:
            
            >>> from Bio import SeqIO
            >>> record = SeqIO.read("GenBank/cor6_6.gb", "genbank")
            Traceback (most recent call last):
                ...
            ValueError: More than one record found in handle
            
            If however you want the first record from a file containing
            multiple records this function would raise an exception (as
            shown in the example above).  Instead use:
            
            >>> from Bio import SeqIO
            >>> record = next(SeqIO.parse("GenBank/cor6_6.gb", "genbank"))
            >>> print("First record's ID %s" % record.id)
            First record's ID X55053.1
            
            Use the Bio.SeqIO.parse(handle, format) function if you want
            to read multiple records from the handle.
        
        to_dict(sequences, key_function=None)
            Turn a sequence iterator or list into a dictionary.
            
            Arguments:
             - sequences  - An iterator that returns SeqRecord objects,
               or simply a list of SeqRecord objects.
             - key_function - Optional callback function which when given a
               SeqRecord should return a unique key for the dictionary.
            
            e.g. key_function = lambda rec : rec.name
            or,  key_function = lambda rec : rec.description.split()[0]
            
            If key_function is omitted then record.id is used, on the assumption
            that the records objects returned are SeqRecords with a unique id.
            
            If there are duplicate keys, an error is raised.
            
            Since Python 3.7, the default dict class maintains key order, meaning
            this dictionary will reflect the order of records given to it. For
            CPython and PyPy, this was already implemented for Python 3.6, so
            effectively you can always assume the record order is preserved.
            
            Example usage, defaulting to using the record.id as key:
            
            >>> from Bio import SeqIO
            >>> filename = "GenBank/cor6_6.gb"
            >>> format = "genbank"
            >>> id_dict = SeqIO.to_dict(SeqIO.parse(filename, format))
            >>> print(list(id_dict))
            ['X55053.1', 'X62281.1', 'M81224.1', 'AJ237582.1', 'L31939.1', 'AF297471.1']
            >>> print(id_dict["L31939.1"].description)
            Brassica rapa (clone bif72) kin mRNA, complete cds
            
            A more complex example, using the key_function argument in order to
            use a sequence checksum as the dictionary key:
            
            >>> from Bio import SeqIO
            >>> from Bio.SeqUtils.CheckSum import seguid
            >>> filename = "GenBank/cor6_6.gb"
            >>> format = "genbank"
            >>> seguid_dict = SeqIO.to_dict(SeqIO.parse(filename, format),
            ...               key_function = lambda rec : seguid(rec.seq))
            >>> for key, record in sorted(seguid_dict.items()):
            ...     print("%s %s" % (key, record.id))
            /wQvmrl87QWcm9llO4/efg23Vgg AJ237582.1
            BUg6YxXSKWEcFFH0L08JzaLGhQs L31939.1
            SabZaA4V2eLE9/2Fm5FnyYy07J4 X55053.1
            TtWsXo45S3ZclIBy4X/WJc39+CY M81224.1
            l7gjJFE6W/S1jJn5+1ASrUKW/FA X62281.1
            uVEYeAQSV5EDQOnFoeMmVea+Oow AF297471.1
            
            This approach is not suitable for very large sets of sequences, as all
            the SeqRecord objects are held in memory. Instead, consider using the
            Bio.SeqIO.index() function (if it supports your particular file format).
            
            This dictionary will reflect the order of records given to it.
        
        write(sequences, handle, format)
            Write complete set of sequences to a file.
            
            Arguments:
             - sequences - A list (or iterator) of SeqRecord objects, or a single
               SeqRecord.
             - handle    - File handle object to write to, or filename as string.
             - format    - lower case string describing the file format to write.
            
            Note if providing a file handle, your code should close the handle
            after calling this function (to ensure the data gets flushed to disk).
            
            Returns the number of records written (as an integer).
    
    FILE
        /home/student/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/__init__.py
    
    



```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```


```python
#we can also use this function on similar files. 

for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
indentfiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
#identifiers will show all the identifiers for every seq record parse chosen. 
 
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    Z78533.1



```python
print(first_record.id)
```

    Z78533.1



```python
print(first_record.description)
```

    C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    Z78532.1



```python
print(second_record.description)
```

    C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("the last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    the last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740



```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78531.1
    Name: Z78531
    Description: C.fasciculatum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78531']
    /sequence_version=1
    /gi=2765656
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium fasciculatum
    /organism=Cypripedium fasciculatum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78531'], 'sequence_version': 1, 'gi': '2765656', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium fasciculatum', 'organism': 'Cypripedium fasciculatum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78531'], 1, '2765656', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium fasciculatum', 'Cypripedium fasciculatum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
 #we can find out the name of the source of the data
print(first_record.annotations["source"])
```

    Cypripedium fasciculatum



```python
#we can also find out the name of the organism
print(first_record.annotations["organism"])
```

    Cypripedium fasciculatum



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum', 'Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'Z78533.1'




```python
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='Z78533', description='C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "desired new description"
```


```python
print(first_record.format("genbank"[:200]))
```

    LOCUS       Z78533                   740 bp    DNA     linear   PLN 30-NOV-2006
    DEFINITION  new_id desired new description.
    ACCESSION   new_id
    VERSION     new_id  GI:2765658
    KEYWORDS    5.8S ribosomal RNA; 5.8S rRNA gene; internal transcribed spacer;
                ITS1; ITS2.
    SOURCE      Cypripedium irapeanum
      ORGANISM  Cypripedium irapeanum
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae;
                Cypripedioideae; Cypripedium.
    REFERENCE   1
      AUTHORS   Cox,A.V., Pridgeon,A.M., Albert,V.A. and Chase,M.W.
      TITLE     Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae):
                nuclear rDNA ITS sequences
      JOURNAL   Unpublished
    REFERENCE   2  (bases 1 to 740)
      AUTHORS   Cox,A.V.
      TITLE     Direct Submission
      JOURNAL   Submitted (19-AUG-1996) Cox A.V., Royal Botanic Gardens, Kew,
                Richmond, Surrey TW9 3AB, UK
    FEATURES             Location/Qualifiers
         source          1..740
                         /organism="Cypripedium irapeanum"
                         /mol_type="genomic DNA"
                         /db_xref="taxon:49711"
         misc_feature    1..380
                         /note="internal transcribed spacer 1"
         gene            381..550
                         /gene="5.8S rRNA"
         rRNA            381..550
                         /gene="5.8S rRNA"
                         /product="5.8S ribosomal RNA"
         misc_feature    551..740
                         /note="internal transcribed spacer 2"
    ORIGIN
            1 cgtaacaagg tttccgtagg tgaacctgcg gaaggatcat tgatgagacc gtggaataaa
           61 cgatcgagtg aatccggagg accggtgtac tcagctcacc gggggcattg ctcccgtggt
          121 gaccctgatt tgttgttggg ccgcctcggg agcgtccatg gcgggtttga acctctagcc
          181 cggcgcagtt tgggcgccaa gccatatgaa agcatcaccg gcgaatggca ttgtcttccc
          241 caaaacccgg agcggcggcg tgctgtcgcg tgcccaatga attttgatga ctctcgcaaa
          301 cgggaatctt ggctctttgc atcggatgga aggacgcagc gaaatgcgat aagtggtgtg
          361 aattgcaaga tcccgtgaac catcgagtct tttgaacgca agttgcgccc gaggccatca
          421 ggctaagggc acgcctgctt gggcgtcgcg cttcgtctct ctcctgccaa tgcttgcccg
          481 gcatacagcc aggccggcgt ggtgcggatg tgaaagattg gccccttgtg cctaggtgcg
          541 gcgggtccaa gagctggtgt tttgatggcc cggaacccgg caagaggtgg acggatgctg
          601 gcagcagctg ccgtgcgaat cccccatgtt gtcgtgcttg tcggacaggc aggagaaccc
          661 ttccgaaccc caatggaggg cggttgaccg ccattcggat gtgaccccag gtcaggcggg
          721 ggcacccgct gagtttacgc
    //
    



```python
from  Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
                 "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
                 "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
                 "SSAC",
                ), 
id = "gi|14150838|gb|AAK54658.1|AF376133_1",
description = "chalcone synthase [Cucumis sativus]",
        
```


```python
rec2 = SeqRecord(
Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
   "DMVVVEIPKLGKEAAVKAIKEWGQ", 
   ),
id = "gi |13919613|gb| AAK33142.1|", 
    description = "chalcone synthase [Fragaria vesca subsp. bracteata]", )

```


```python
rec3 = SeqRecord(
Seq("MTHSAIUHDAOSIJAOIHSCBASLHPPPDJALSKJALISNCLIS"
   "EKSHAKFHEINAISHFOAINFEFPOEJFPOEHFWIEGFOWEIHWO"
   "KASFHWEUHOWINEFNWUORFOLEIQEHOIHIHOQHIWPQPQPQP"
   "HHIHFQOEHWOENFWOEUGORUGBWOUDBBCBCNSSIWPPPQJEK"
   "QQWWEDRTYUUFDIJANKAAPPBCHSKWIRHNSSWJJNLQPEJDM"
   "QQBJHBFBWPEIBLASKDNIDNQDNOIF", 
   ), 
    id = "gi|13925890|gb|AAK49457.1|",
    description= "chalcone synthase [Nicotiana tabacum]",)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa","fasta" )
```


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-102-d88394100fac> in <module>
          1 from Bio import SeqIO
    ----> 2 SeqIO.write(my_records, "my_example.faa","fasta" )
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/__init__.py in write(sequences, handle, format)
        516         with as_handle(handle, "w") as fp:
        517             for record in sequences:
    --> 518                 fp.write(format_function(record))
        519                 count += 1
        520         return count


    ~/anaconda3/lib/python3.7/site-packages/Bio/SeqIO/FastaIO.py in as_fasta(record)
        414     method and by the SeqIO.write(..., ..., "fasta") function.
        415     """
    --> 416     id = _clean(record.id)
        417     description = _clean(record.description)
        418     if description and description.split(None, 1)[0] == id:


    AttributeError: 'tuple' object has no attribute 'id'



```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count= SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(records)
```

    [SeqRecord(seq=Seq('GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGT...ACG'), id='rc_Z78533.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAAT...ACG'), id='rc_Z78532.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCC...ACG'), id='rc_Z78531.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGG...ACG'), id='rc_Z78530.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAAC...CGT'), id='rc_Z78529.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGT...ACG'), id='rc_Z78527.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAAT...ACG'), id='rc_Z78526.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATC...ACA'), id='rc_Z78525.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAAC...ACG'), id='rc_Z78524.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCA...ACG'), id='rc_Z78523.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCC...ACG'), id='rc_Z78522.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAG...TAC'), id='rc_Z78521.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACC...ACG'), id='rc_Z78520.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAA...TAT'), id='rc_Z78519.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATC...ACG'), id='rc_Z78518.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAAT...ACG'), id='rc_Z78517.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCA...ACG'), id='rc_Z78516.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAAT...ACG'), id='rc_Z78515.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCC...ACG'), id='rc_Z78514.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGG...ACG'), id='rc_Z78513.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGT...ACG'), id='rc_Z78512.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGA...ACG'), id='rc_Z78511.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...TAG'), id='rc_Z78510.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACC...ACG'), id='rc_Z78509.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGG...ACG'), id='rc_Z78508.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGG...ACG'), id='rc_Z78507.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCC...ACG'), id='rc_Z78506.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...ACG'), id='rc_Z78505.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCA...ACG'), id='rc_Z78504.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...ACG'), id='rc_Z78503.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGT...ACG'), id='rc_Z78502.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGT...ACG'), id='rc_Z78501.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTC...ACG'), id='rc_Z78500.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGT...ACG'), id='rc_Z78499.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAAC...ACG'), id='rc_Z78498.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_Z78497.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCA...ACG'), id='rc_Z78496.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGA...ACG'), id='rc_Z78495.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATT...ACG'), id='rc_Z78494.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATG...ACG'), id='rc_Z78493.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAAT...ACG'), id='rc_Z78492.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAA...ACG'), id='rc_Z78491.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCA...ACG'), id='rc_Z78490.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_Z78489.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGG...CAG'), id='rc_Z78488.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATC...ACG'), id='rc_Z78487.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGAT...ACG'), id='rc_Z78486.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACT...CAG'), id='rc_Z78485.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCA...ACG'), id='rc_Z78484.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_Z78483.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...AGA'), id='rc_Z78482.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGAT...ACG'), id='rc_Z78481.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGAT...ACG'), id='rc_Z78480.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTG...ACG'), id='rc_Z78479.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAA...ACG'), id='rc_Z78478.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTC...ACG'), id='rc_Z78477.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGG...ACG'), id='rc_Z78476.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCC...ACG'), id='rc_Z78475.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAAC...ACG'), id='rc_Z78474.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAAT...ACG'), id='rc_Z78473.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_Z78472.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_Z78471.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTAT...ACG'), id='rc_Z78470.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTA...ACG'), id='rc_Z78469.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAA...ACG'), id='rc_Z78468.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATC...ACG'), id='rc_Z78467.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGG...ACG'), id='rc_Z78466.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCA...ACG'), id='rc_Z78465.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCA...ACG'), id='rc_Z78464.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTC...ACG'), id='rc_Z78463.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAAT...ACG'), id='rc_Z78462.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCA...ACG'), id='rc_Z78461.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCA...ACG'), id='rc_Z78460.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAAC...ACG'), id='rc_Z78459.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAAC...ACG'), id='rc_Z78458.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGA...ACG'), id='rc_Z78457.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAA...ACG'), id='rc_Z78456.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCA...ACG'), id='rc_Z78455.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGG...ACG'), id='rc_Z78454.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_Z78453.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATC...ACG'), id='rc_Z78452.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCA...ACG'), id='rc_Z78451.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCA...TCC'), id='rc_Z78450.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_Z78449.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTG...ACG'), id='rc_Z78448.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_Z78447.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGG...ACG'), id='rc_Z78446.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCA...ACG'), id='rc_Z78445.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCT...ACG'), id='rc_Z78444.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGG...ACG'), id='rc_Z78443.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATC...TAC'), id='rc_Z78442.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCA...TCC'), id='rc_Z78441.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_Z78440.1', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATA...ATG'), id='rc_Z78439.1', name='<unknown name>', description='reverse complement', dbxrefs=[])]



```python
len(records)
```




    94




```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
    if len(rec) < 700
]
```


```python
len(records)
```




    18




```python
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
    if len(rec) < 700 
)

SeqIO.write(records, "rev_complement.fasta", "fasta")
```




    18




```python

```
# Sequence Alignment 

In this module, we identified and aligned sequecnes. We were introduced to the parse function as well. 

#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05371_seed.sth
```


```python
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print("Alignment lnegth %i" % alignment.get_alignment_length())
```

    Alignment lnegth 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
for record in alignment :
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
align1 = MultipleSeqAlignment(
...     [
...      SeqRecord(Seq("ACTGCTAGCTAG"), id = "Alpha"),
...      SeqRecord(Seq("ACT-CTAGCTAG"), id = "Beta"),
...      SeqRecord(Seq("ACTGCTAGDTAG"), id = "Gamma"), 
...     ]
... )
>>> align2 = MultipleSeqAlignment(
...     [ 
...          SeqRecord(Seq("GTCAGC-AG"), id = "Delta"), 
...          SeqRecord(Seq("GACAGCTAG"), id = "Epsilon"), 
...          SeqRecord(Seq("GTCAGCTAG"), id = "Zeta"),
...     ]
... )
>>> align3 = MultipleSeqAlignment(
...      [
...          SeqRecord(Seq("ACTAGTACAGCTG"), id = "Eta", )
...      ]
... )
```


```python
my_alignments = [align1, align2, align3]
```


```python
#this will tell us what is in our alignments 

my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7fb7db5d9350>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7fb7db5d9450>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 13) at 7fb7db5d94d0>]




```python
#we can also priint out alignments
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7fb7db5d9350>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7fb7db5d9450>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 13) at 7fb7db5d94d0>]



```python
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 1 rows and 13 columns
    ACTAGTACAGCTG Eta
    



```python
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
last_align = alignments[-1]
```


```python
#we can print our last group of alignments 
print(last_align)
```

    Alignment with 1 rows and 13 columns
    ACTAGTACAGCTG Eta



```python
#we can also print our first group of alignments
first_align = alignments [0]
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
from Bio import AlignIO
```


```python
count = AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth.txt", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" %i

```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print("Number of rows: %i " % len(alignment))
```

    Number of rows: 7 



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
print(alignment[2,6])
```

    T



```python
print(alignment[2].seq[6])
```

    T



```python
print(alignment[:, 6])
```

    TTT---T



```python
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited.sort()
```


```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Align import MultipleSeqAlignment
```


```python
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id = "seq1"),
        SeqRecord(Seq("AAT-CTA"), id = "seq2"),
        SeqRecord(Seq("CCTACT-"), id = "seq3"),
        SeqRecord(Seq("TCTCCTC"), id = "seq4"),
    ]
)
```


```python
substitutions = alignment.substitutions
```


```python
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
m = substituions.select("ACTG")
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
help(ClustalwCommandline)
```

    Help on class ClustalwCommandline in module Bio.Align.Applications._Clustalw:
    
    class ClustalwCommandline(Bio.Application.AbstractCommandline)
     |  ClustalwCommandline(cmd='clustalw', **kwargs)
     |  
     |  Command line wrapper for clustalw (version one or two).
     |  
     |  http://www.clustal.org/
     |  
     |  Notes
     |  -----
     |  Last checked against versions: 1.83 and 2.1
     |  
     |  References
     |  ----------
     |  Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA,
     |  McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD,
     |  Gibson TJ, Higgins DG. (2007). Clustal W and Clustal X version 2.0.
     |  Bioinformatics, 23, 2947-2948.
     |  
     |  Examples
     |  --------
     |  >>> from Bio.Align.Applications import ClustalwCommandline
     |  >>> in_file = "unaligned.fasta"
     |  >>> clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
     |  >>> print(clustalw_cline)
     |  clustalw2 -infile=unaligned.fasta
     |  
     |  You would typically run the command line with clustalw_cline() or via
     |  the Python subprocess module, as described in the Biopython tutorial.
     |  
     |  Method resolution order:
     |      ClustalwCommandline
     |      Bio.Application.AbstractCommandline
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self, cmd='clustalw', **kwargs)
     |      Initialize the class.
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from Bio.Application.AbstractCommandline:
     |  
     |  __call__(self, stdin=None, stdout=True, stderr=True, cwd=None, env=None)
     |      Execute command, wait for it to finish, return (stdout, stderr).
     |      
     |      Runs the command line tool and waits for it to finish. If it returns
     |      a non-zero error level, an exception is raised. Otherwise two strings
     |      are returned containing stdout and stderr.
     |      
     |      The optional stdin argument should be a string of data which will be
     |      passed to the tool as standard input.
     |      
     |      The optional stdout and stderr argument may be filenames (string),
     |      but otherwise are treated as a booleans, and control if the output
     |      should be captured as strings (True, default), or ignored by sending
     |      it to /dev/null to avoid wasting memory (False). If sent to a file
     |      or ignored, then empty string(s) are returned.
     |      
     |      The optional cwd argument is a string giving the working directory
     |      to run the command from. See Python's subprocess module documentation
     |      for more details.
     |      
     |      The optional env argument is a dictionary setting the environment
     |      variables to be used in the new process. By default the current
     |      process' environment variables are used. See Python's subprocess
     |      module documentation for more details.
     |      
     |      Default example usage::
     |      
     |          from Bio.Emboss.Applications import WaterCommandline
     |          water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
     |                                       stdout=True, auto=True,
     |                                       asequence="a.fasta", bsequence="b.fasta")
     |          print("About to run: %s" % water_cmd)
     |          std_output, err_output = water_cmd()
     |      
     |      This functionality is similar to subprocess.check_output(). In general
     |      if you require more control over running the command, use subprocess
     |      directly.
     |      
     |      When the program called returns a non-zero error level, a custom
     |      ApplicationError exception is raised. This includes any stdout and
     |      stderr strings captured as attributes of the exception object, since
     |      they may be useful for diagnosing what went wrong.
     |  
     |  __repr__(self)
     |      Return a representation of the command line object for debugging.
     |      
     |      e.g.
     |      
     |      >>> from Bio.Emboss.Applications import WaterCommandline
     |      >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
     |      >>> cline.asequence = "asis:ACCCGGGCGCGGT"
     |      >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
     |      >>> cline.outfile = "temp_water.txt"
     |      >>> print(cline)
     |      water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
     |      >>> cline
     |      WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5)
     |  
     |  __setattr__(self, name, value)
     |      Set attribute name to value (PRIVATE).
     |      
     |      This code implements a workaround for a user interface issue.
     |      Without this __setattr__ attribute-based assignment of parameters
     |      will silently accept invalid parameters, leading to known instances
     |      of the user assuming that parameters for the application are set,
     |      when they are not.
     |      
     |      >>> from Bio.Emboss.Applications import WaterCommandline
     |      >>> cline = WaterCommandline(gapopen=10, gapextend=0.5, stdout=True)
     |      >>> cline.asequence = "a.fasta"
     |      >>> cline.bsequence = "b.fasta"
     |      >>> cline.csequence = "c.fasta"
     |      Traceback (most recent call last):
     |      ...
     |      ValueError: Option name csequence was not found.
     |      >>> print(cline)
     |      water -stdout -asequence=a.fasta -bsequence=b.fasta -gapopen=10 -gapextend=0.5
     |      
     |      This workaround uses a whitelist of object attributes, and sets the
     |      object attribute list as normal, for these.  Other attributes are
     |      assumed to be parameters, and passed to the self.set_parameter method
     |      for validation and assignment.
     |  
     |  __str__(self)
     |      Make the commandline string with the currently set options.
     |      
     |      e.g.
     |      
     |      >>> from Bio.Emboss.Applications import WaterCommandline
     |      >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
     |      >>> cline.asequence = "asis:ACCCGGGCGCGGT"
     |      >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
     |      >>> cline.outfile = "temp_water.txt"
     |      >>> print(cline)
     |      water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
     |      >>> str(cline)
     |      'water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5'
     |  
     |  set_parameter(self, name, value=None)
     |      Set a commandline option for a program (OBSOLETE).
     |      
     |      Every parameter is available via a property and as a named
     |      keyword when creating the instance. Using either of these is
     |      preferred to this legacy set_parameter method which is now
     |      OBSOLETE, and likely to be DEPRECATED and later REMOVED in
     |      future releases.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Bio.Application.AbstractCommandline:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Bio.Application.AbstractCommandline:
     |  
     |  parameters = None
    



```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.aln
```


```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln.txt", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
tree = Phylo.read("opuntia.dnd.txt", "newick")
```


```python
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    



```python

```
# Blast

In this module, we learned how to run a BLAST from importing a fasta file. 

from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "Dr. Josh.vandenbrink@gmail.com"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
#https://github.com/biopython/biopython/blob/master/Doc/examples/m_cold.fasta
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta", format = "fasta")
```


```python
print(record)
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
with open ("m_cold.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```


```python

```


```python

```


```python

```
# Challenge 1

In this challenge, we ran a BLAST on a cancer gene of our choice. 

from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "sel035@email.latech.edu"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "64002")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("sequence.fasta", format = "fasta")
```


```python
print(record)
```

    ID: NC_000002.12:192749845-192776899
    Name: NC_000002.12:192749845-192776899
    Description: NC_000002.12:192749845-192776899 Homo sapiens chromosome 2, GRCh38.p14 Primary Assembly
    Number of features: 0
    Seq('AAGGCACTCTGGCACCCAGTTTTGGAACTGCAGTTTTAAAAGTCATAAATTGAA...GAA')



```python
result_handle = NCBIWWW.qblast("blastn", "8332116", record.seq)
```


```python
with open ("sequence.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("sequence.fasta")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```


```python

```
# Open CV

In this module, we learned how to open pictures, flip them, and change their color. 

import numpy as np 
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
img = cv2.imread("mushroom.jpeg")
```


```python
type(img)
```




    numpy.ndarray




```python

```


```python
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
# we can use this to print the image

plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f4ab4b7d310>




![png](output_7_1.png)



```python
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7f4ab22ab5d0>




![png](output_9_1.png)



```python
img_gray = cv2.imread("mushroom.jpeg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (398, 600)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7f4ab1973ad0>




![png](output_11_1.png)



```python
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f4ab18e1550>




![png](output_12_1.png)



```python
fix_img.shape
```




    (398, 600, 3)




```python
new_img = cv2.resize(fix_img, (1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f4ab18ca550>




![png](output_14_1.png)



```python
new_img.shape
```




    (400, 1000, 3)




```python
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f4ab18346d0>




![png](output_17_1.png)



```python
new_img.shape
```




    (400, 1000, 3)




```python
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7f4ab1795910>




![png](output_19_1.png)



```python
#we can also flip the image on a different axis/different way

flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7f4ab16f1790>




![png](output_20_1.png)



```python
type(fix_img)
```




    numpy.ndarray




```python
cv2.imwrite('Mushroom_fixed_image.jpg', fix_img)
```




    True




```python

```
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("mushroom.jpeg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f3ab4153650>




![png](output_2_1.png)



```python
#we can fix the color scale
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f3aafc595d0>




![png](output_4_1.png)



```python

```


```python

```


```python

```


```python
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f3aaccc8fd0>




![png](output_9_1.png)



```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f3aaccbe790>




![png](output_11_1.png)



```python
img1 = cv2.imread('do-not-copy-stamp.jpeg')
img2 = cv2.imread('mushroom.jpeg')
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f3aacc2ae10>




![png](output_13_1.png)



```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f3aacba14d0>




![png](output_15_1.png)



```python
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200,1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


    ---------------------------------------------------------------------------

    error                                     Traceback (most recent call last)

    <ipython-input-28-fbcc42a1e8a2> in <module>
    ----> 1 blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
    

    error: OpenCV(4.9.0) /io/opencv/modules/core/src/arithm.cpp:647: error: (-209:Sizes of input arguments do not match) The operation is neither 'array op array' (where arrays have the same size and the same number of channels), nor 'array op scalar', nor 'scalar op array' in function 'arithm_op'




```python
plt.imshow(blended)
```


```python
alpha = 0.8
beta = 0.2

blended = cv2.addWeighted(img2, alpha, img1, beta, 0)
plt.imshow(blended)
```


```python
img1 = cv2.imread('do-not-copy-stamp.jpeg')
img2= cv2.imread('mushroom.jpeg')



img1 = cv2.resize(img1, (600,600))
```


```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
large_img = img2
small_img = img1

x_offset= 0
y_offset= 0

x_end= x_offset + small_img.shape[1]
y_end= y_offset + small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```


```python

```


```python

```


```python

```


```python

```
#https://github.com/worklifesg/Python-for-Computer-Vision-with-OpenCV-and-Deep-Learning
```


```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline

```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f848950fd50>




![png](output_3_1.png)



```python
img = cv2.imread('rainbow.jpg', 0)
```


```python
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f8489438510>




![png](output_5_1.png)



```python
ret1, thresh1 = cv2.threshold(img, 127,255, cv2.THRESH_BINARY)
```


```python
ret1

```




    127.0




```python
plt.imshow(thresh1, cmap="gray")
```




    <matplotlib.image.AxesImage at 0x7f848839d5d0>




![png](output_8_1.png)



```python
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap ="gray")
```




    <matplotlib.image.AxesImage at 0x7f8488181a50>




![png](output_9_1.png)



```python

```


```python
img3 = cv2.imread('rainbow.jpg')
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap ="gray")
```




    <matplotlib.image.AxesImage at 0x7f84880e65d0>




![png](output_11_1.png)



```python
img_r = cv2.imread('crossword.jpg', 0)
```


```python
def show_pic(img):
    fig= plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```


![png](output_14_0.png)



```python
ret, th1 =cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_15_0.png)



```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_16_0.png)



```python
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY,11,8)
```


```python
show_pic(th2)
```


![png](output_18_0.png)



```python
#we can blend this image into new thresholds

blended = cv2.addWeighted(src1= th1, alpha =0.6,
                        src2 = th2, beta = 0.4, gamma=0)
```


```python
show_pic(blended)
```


![png](output_20_0.png)



```python
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11,8)

blended = cv2.addWeighted(src1= th1, alpha =0.6,
                        src2 = th2, beta = 0.4, gamma=0)
show_pic(blended)
```


![png](output_21_0.png)



```python

```

# Aspect Detection: Corner detection

In this module, we learned about the Shi-Tomasai Corner detection.
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
flat_chess = cv2.imread('green_chessboard.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fcab7b8f610>




![png](output_1_1.png)



```python
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fcab72ead90>




![png](output_2_1.png)



```python
real_chess = cv2.imread("chessboard.jpeg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fcab725d210>




![png](output_4_1.png)



```python
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fcab4005110>




![png](output_5_1.png)



```python
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize =3, k =0.04)
dst = cv2.dilate(dst, None)
```


```python
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fcaaf7ed510>




![png](output_7_1.png)



```python
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize =3, k =0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fcaaf745f90>




![png](output_8_1.png)



```python
#Shi-Tomasai Corner Detection

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y), 3, (255,0,0), -1)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fcaaf739310>




![png](output_10_1.png)



```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (255,0,0), -1)
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fcaaf691810>




![png](output_11_1.png)



```python

```
# Aspect Detection: Edge Detection

import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("mushroom.jpeg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fa77fd67a50>




![png](output_3_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77edbca50>




![png](output_4_1.png)



```python
med_value = np.median(img)
med_value
```




    95.0




```python
lower = int(max(0, 0.7*med_value))
upper = int(min(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c538250>




![png](output_6_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c49ead0>




![png](output_7_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1=lower,
                 threshold2=upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c39d390>




![png](output_8_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1=lower,
                 threshold2=upper +50)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c387090>




![png](output_9_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1=lower,
                 threshold2=upper +100)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c2e6dd0>




![png](output_10_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1=lower,
                 threshold2=upper +60)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c251ad0>




![png](output_11_1.png)



```python
blurred_img = cv2.blur(img, ksize = (8,8))

edges = cv2.Canny(image=blurred_img,
                 threshold1=lower,
                 threshold2=upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fa77c2407d0>




![png](output_12_1.png)



```python

```
# Feauture Detection: Feature Matches

In this module, we looked at the similarities between two pictures and found their matches. 

import cv2
import numpy as np 
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display (img, cmap = 'gray'):
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
cheerios = cv2.imread("cheerios.jpeg")
display(cheerios)
```


![png](output_2_0.png)



```python
cereals = cv2.imread('cereals.jpeg',0)
display(cereals)
```


![png](output_3_0.png)



```python
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(cheerios, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1,des2)
```


```python
matches = sorted(matches, key = lambda x:x.distance)
```


```python
cheerios_matches = cv2.drawMatches(cheerios, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
display(cheerios_matches)
```


![png](output_8_0.png)



```python
sift = cv2.SIFT_create()
```


```python
kp1,des1 = sift.detectAndCompute(cheerios, None)
kp2,des2 = sift.detectAndCompute(cereals, None)
```


```python
bf = cv2.BFMatcher()
matches= bf.knnMatch(des1, des2, k=2)
```


```python
good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('Length of total matches:', len(matches))
print("length of good matches", len(good))
```

    Length of total matches: 631
    length of good matches 35



```python
sift_matches = cv2.drawMatchesKnn(cheerios, kp1, cereals, kp2, good, None, flags = 2 )
display(sift_matches)
```


![png](output_14_0.png)



```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(cheerios, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees=5)
search_params = dict(checks=50)

```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
flann_matches = cv2.drawMatchesKnn(cheerios, kp1, cereals, kp2, good, None, flags =0)
display(flann_matches)
```


![png](output_18_0.png)



```python
sift = cv2.SIFT_create()
kp1, des1 = sift.detectAndCompute(cheerios, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees =5)
search_param = dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_param)

matches = flann.knnMatch(des1, des2, k=2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0,255,0),
                   singlePointcolor = (255,0,0), 
                   matchesMask = matchesMask,
                   flags =0)
```


```python
flann_matches = cv2.drawMatchesKnn(cheerios, kp1, cereals, kp2, matches, None, **draw_params)

display(flann_matches)
```


    ---------------------------------------------------------------------------

    error                                     Traceback (most recent call last)

    <ipython-input-36-88bd83582a92> in <module>
    ----> 1 flann_matches = cv2.drawMatchesKnn(cheerios, kp1, cereals, kp2, matches, None, **draw_params)
          2 
          3 display(flann_matches)


    error: OpenCV(4.9.0) :-1: error: (-5:Bad argument) in function 'drawMatchesKnn'
    > Overload resolution failed:
    >  - 'singlePointcolor' is an invalid keyword argument for drawMatchesKnn()
    >  - 'singlePointcolor' is an invalid keyword argument for drawMatchesKnn()




```python

```
# Feature Detection: Object Detection

In this module, we looked at two images to build a heat map and specifically pick out the object we trained it to. 

import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline 
```


```python
full = cv2.imread('testing_sunflower.jpeg')
```


```python
full =cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7f35b87b5f10>




![png](output_3_1.png)



```python
full2 = cv2.imread('training_sunflower.jpeg')
```


```python
full2 =cv2.cvtColor(full2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full2)
```




    <matplotlib.image.AxesImage at 0x7f35b867aa10>




![png](output_6_1.png)



```python
test = cv2.imread('testing_sunflower.jpeg')
```


```python
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7f35b4bff410>




![png](output_9_1.png)



```python
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (900, 891, 3)
    Training image shape: (900, 891, 3)



```python
methods =['cv2.TM_CCOEFF','cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
for m in methods:
    test_copy = test.copy()
    method = eval(m)
    
    res= cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc, = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
        
        height, width, channels = full.shape
        bottom_right = (top_left[0] + width, top_left[1] +height)
        
        cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0), 10)
        
        plt.subplot(121)
        plt.imshow(res)
        plt.title("heatmap of template matching")
        plt.subplot(122)
        plt.imshow(test_copy)
        plt.title('Detection of template')
        
        plt.suptitle(m)
        
        plt.show()
        print('\n')
        print('\n')
```


![png](output_12_0.png)


    
    
    
    



![png](output_12_2.png)


    
    
    
    



![png](output_12_4.png)


    
    
    
    



![png](output_12_6.png)


    
    
    
    



```python

```


```python

```
