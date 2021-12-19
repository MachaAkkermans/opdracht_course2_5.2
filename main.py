import re
from translation import code


class Sequentie:

    def __init__(self, sequentie):
        self.__sequentie = sequentie

    def setsequentie(self, sequentie):
        self.__sequentie = sequentie

    def getlengte(self):
        self.__lengte = len(self.__sequentie)
        return self.__lengte

    def getsequentie(self):
        return self.__sequentie


class Rna(Sequentie):

    def __init__(self, sequentie):
        self.setrna(sequentie)
        # self.gettranslatie()
        super().__init__(self.__sequentie)

    def setrna(self, sequentie):
        sequentie = sequentie.replace(" ", "")
        x = bool(re.search("^[autcgn]+$", sequentie))
        if x == True:
            self.__sequentie = sequentie
        else:
            print("Dit is geen RNA sequentie")
            self.__sequentie = ""

    def getrna(self):
        return self.__sequentie

    def gettranslatie(self):
        self.__eiwit = ""
        self.__sequentie = self.__sequentie.lower()
        print(self.__sequentie)
        codon = ""
        for i in self.__sequentie:
            if len(codon) != 3:
                codon += i
            else:
                if codon in code:
                    # lijst_met_eiwitten.append(code[i])
                    self.__eiwit += code[codon]
                codon = ""
        return self.__eiwit


class Dna(Sequentie):

    def __init__(self, sequentie):
        self.setdna(sequentie)
        super().__init__(self.__sequentie)

    def setdna(self, sequentie):
        x = bool(re.search("^[ATGCN]+$", sequentie))
        if x == True:
            self.__sequentie = sequentie
        else:
            print("dit is geen DNA sequentie")
            self.__sequentie = ""

    def getdna(self):
        return self.__sequentie

    def gettranscript(self):
        self.__transcript = ""
        for letter in self.__sequentie:
            if letter == 'A':
                self.__transcript += 'U'
            elif letter == 'T':
                self.__transcript += 'A'
            elif letter == 'G':
                self.__transcript += 'C'
            else:
                self.__transcript += 'G'

        self.__transcript = self.__transcript[::-1]
        return self.__transcript

    def getgc(self):
        try:
            totaal = len(self.__sequentie)
            aantal_g = self.__sequentie.count("G")
            aantal_c = self.__sequentie.count("C")
            gc = int(round(((aantal_c + aantal_g) / totaal * 100), 2))
            return gc
        except ZeroDivisionError:
            return ""


def readfile1(filename):
    file = open(filename, "r")
    lijst_met_dna = []
    inhoud = file.read()
    mix = re.split(">", inhoud, re.MULTILINE)
    # print(mix)
    del mix[0]
    # print(mix)
    for fasta in mix:
        header_met_sequentie = []
        header, sequentie = fasta.split("\n", 1)
        sequentie = sequentie.replace("\n", "")
        header_met_sequentie.append(header)
        header_met_sequentie.append(sequentie)
        lijst_met_dna.append(header_met_sequentie)
    file.close()
    # print(lijst_met_dna)
    return lijst_met_dna


def rnaaanvragen(lijst_met_rna):
    sequentie = "at gggatcaaca acaatgtccc ctccatcttttcccgtcgtc ctcctgctcc " \
                "tcctcctcgc caccatagcc gcagccgccg gaagcaacatggatgaggag gtggtggacg" \
                "acctccagta tcttattgac aactccgacg acatccccaccaacgatccc gacgggtggc " \
                "ctgagggaga ctacgacgac gacgaccttc tcttccaagatcaggaccag gacctcacag " \
                "gccaccagcc ggagatcgac gagacccacg tcgtggtcctcgccgccgca aacttttcct " \
                "ccttcctcgc ctccagccac catgttatgg ttgagttctacgcaccttgg tgtggccact" \
                "gccaggagct cgccccggat tacgccgccg ccgccgcgcatctcgccgct caccaccacc" \
                " aggcccacct cgcccttgcc aaggtcgacg ccaccgagga"

    r = Rna(sequentie)
    eiwit = r.gettranslatie()
    print(eiwit)


def dnaaanvragen(lijst_met_dna):
    gc_hoogste = 0
    lijstdna = []
    for fasta in lijst_met_dna:
        d = Dna(fasta[1])
        if d.getdna() != "":
            lijstdna.append(d.getdna())
            gc = d.getgc()
            if gc >= gc_hoogste:
                sequentie_langste = fasta[1]
                lengte = (d.getlengte())
                rna = (d.gettranscript())
                gc_hoogste = gc
                gene_naam_locatie = fasta[0].find("gene:")
                gene_naam_locatie_stop = fasta[0][gene_naam_locatie:] \
                    .find(" ")
                gene_naam = fasta[0][gene_naam_locatie:(
                        gene_naam_locatie_stop + gene_naam_locatie)]
    print("gen naam = ", gene_naam)
    print("sequentie = ", sequentie_langste)
    print("GC% =", gc_hoogste)
    print("lengte = ", lengte)
    print("rna = ", rna)
    print()
    lijstdna = sorted(lijstdna, key=len)
    print(lijstdna)






def main():
    # bestandsnaam = input("Voer je bestandsnaam in")
    bestandsnaam = "felis_catus.fa"
    # bestandsnaam = "new.txt"
    # bestandsnaam = "bestand.txt"
    # bestandsnaam = "test.txt"
    lijst_met_sequentie = readfile1(bestandsnaam)
    dnaaanvragen(lijst_met_sequentie)
    #rnaaanvragen(lijst_met_sequentie)



main()
