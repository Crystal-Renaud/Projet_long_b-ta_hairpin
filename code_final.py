import os
import sys
import pandas as pd
import itertools 

#extraire le nombre d'acide aminé présent dans la protéine
def nb_AA(prot):
    with open(prot, 'r') as f: 
        for i in range(4): #pour lire les 4 premieres lignes
            line = f.readline() #lis les lignes une par une
            list_word = line.split() 
            if len(list_word) == 4: #quand le nb de mot sur une ligne =4
                nb_residu = int(list_word[-1]) #nb retient le dernier mot de la ligne
        return (nb_residu)


#création des listes avec positions, résidus et structures de la protéine
def crea_list(table):
    with open (table, 'r') as f:
        position = [] #liste vide pour les position
        residu = [] #pour les résidus
        structure = [] #pour les structures secondaires
        for line in f:
            pos = line.split()[0] #lis la 1ere colonne
            position.append(pos)
            res = line.split()[4] #lis colonne 5
            residu.append(res)
            struc = line[23:24] #lis colonne 7 et transforme espace en c (coil)
            if struc == ' ':
                struc = "c"
            structure.append(struc)
        return(position, residu, structure)    


#lecture du fichier .hpin 
def read_hpin(hpin):
    with open(hpin, 'r') as file:
        lines = []
        for line in file:
            line = line.rstrip()
            lines.append(line)
        return lines[1:len(lines)]


#Pour récupérer les séquences des différentes listes obtenues après extraction du fichier .hpin
def extraire_si_pos_dans_liste_existe(nom_liste, position):
    if position <= len(nom_liste):
        return nom_liste[position]
    else:
        return ""


# passer le fichiers
print("Entrez un fichier sst obtenu avec Promotif")
prot = input()
if not os.path.exists(prot):
    sys.exit(f"ERREUR: {prot} n'existe pas.")
table = input("Entrez le nom de la protéine traitée et ajouter .txt : ")

print("Entrez le fichier .hpin obtenu avec Promotif")
hpin = input()
if not os.path.exists(hpin):
    sys.exit(f"ERREUR: {hpin} n'existe pas.")

nb = nb_AA(prot)

with open(prot, 'r') as f:
    #lire et stocker toutes les lignes 
    lines = f.readlines()

with open(table, "w") as f:
    for line in lines:
        #le pointeur part du début de la ligne
        f.seek(0)
        #tronque le fichier
        f.truncate()

        #on commence à écrire à la ligne 7 et on s'arrête à nb+7
        f.writelines(lines[7:(nb+7)])
        
position, residu, structure = crea_list(table)

test = read_hpin(hpin)

print( )
print(f"""Les positions associées à la séquence sont: 
{position}""")
print(f"""La séquence de la protéine est: 
{residu}""")
print(f"""Les structures associées à la séquence sont: 
{structure}""")
print(  )
print(  )
print( )

pos_brin_1 = []
seq_brin_1 = []
hairpin = []
pos_brin_2 = []
seq_brin_2 = []
type_hpin = []
type_b_turn = []
rama = []

for j in range(len(test)//6):

    nbLines = j*6

    #extraction position de début et de fin du brin 1
    pos_brin_1.append(test[nbLines][10:13])
    pos_brin_1.append(test[nbLines][17:20])
    #print(f"position début et fin du brin 1 est {pos_brin_1}")

    #extraction séquence du brin 1 du feuillet 
    seq_brin_1.append(test[nbLines+1])
    #print(f"séquence du brin 1 est {seq_brin_1}")

    #extraction séquence du hairpin
    hairpin.append(test[nbLines+2])
    #print(f"hairpin est {hairpin}")

    #extraction position de début et de fin du brin 2 du feuillet
    pos_brin_2.append(test[nbLines][24:27])
    pos_brin_2.append(test[nbLines][31:34])
    #print(f"position début et fin du brin 2 est {pos_brin_2}")
    
    #extraction séquence du brin 2
    seq_brin_2.append(test[nbLines+3])
    #print(f"séquence du brin 2 est {seq_brin_2}")

    #extraction du type de hairpin
    type_hpin.append(test[nbLines][37:43])
    #print(f"le type du hpin est {type_hpin}")

    #extraction du type de beta turn
    for i in (test[nbLines][45:47]):
        if i == "  ":
            " "
        else:
            type_b_turn.append(test[nbLines][45:47])
    #print(type_b_turn)

    #extraction de la zone de Ramachandran du b-turn
    rama.append(test[nbLines+4])
    #print(f"la zone de Ramachandran est {rama}")


#Pour le premier feuillet = bloc a
seq1a = extraire_si_pos_dans_liste_existe(seq_brin_1, 0)
seq1a = seq1a.replace(" ","")
hairpin_a = extraire_si_pos_dans_liste_existe(hairpin, 0)
seq2a = extraire_si_pos_dans_liste_existe(seq_brin_2, 0)
pos_seq1a = extraire_si_pos_dans_liste_existe(pos_brin_1, 0) + ", " + extraire_si_pos_dans_liste_existe(pos_brin_1, 1)
pos_seq2a = extraire_si_pos_dans_liste_existe(pos_brin_2, 0) + ", " + extraire_si_pos_dans_liste_existe(pos_brin_2, 1)
type_hairpin_a = extraire_si_pos_dans_liste_existe(type_hpin, 0)
type_bturn_a = extraire_si_pos_dans_liste_existe(type_b_turn, 0)
zone_rama_a = extraire_si_pos_dans_liste_existe(rama, 0)

print(f"Les positions de début et de fin du brin 1 sont {pos_seq1a}")
print(f"La séquence du brin 1 est {seq1a}")
print(f"Les positions de début et de fin du brin 2 sont {pos_seq2a}")
print(f"La séquence du brin 2 est {seq2a}")
print(f"La séquence du b-hairpin est {hairpin_a}")
print(f"Le type de b-hairpin est {type_hairpin_a}")
print(f"Le type de b-turn est {type_bturn_a}")
print(f"La zone de Ramachandrant est {zone_rama_a}")
print(  )
print( )
print( )

#Vérifier si les séquences du fichier hpin existent dans la protéine
seq_proteine = ''.join(residu)
#pour seq brin 1
if seq1a in seq_proteine:
    print(f"Le brin {seq1a} existe dans la protéine")
else:
    print(f"Le brin {seq1a} n'est pas dans la protéine")

#pour seq brin 2
if seq2a in seq_proteine:
    print(f"Le brin {seq2a} existe dans la protéine")
else:
    print(f"Le brin {seq2a} n'est pas dans la protéine")

#pour seq hairpin
if hairpin_a in seq_proteine:
    print(f"Le brin {hairpin_a} existe dans la protéine")
else:
    print(f"Le brin {hairpin_a} n'est pas dans la protéine")

print()

#Pour le deuxième feuillet = bloc b
seq1b = extraire_si_pos_dans_liste_existe(seq_brin_1, 1)
seq1b = seq1b.replace(" ","")
hairpin_b = extraire_si_pos_dans_liste_existe(hairpin, 1)
seq2b = extraire_si_pos_dans_liste_existe(seq_brin_2, 1)
pos_seq1b = extraire_si_pos_dans_liste_existe(pos_brin_1, 2) + ", " + extraire_si_pos_dans_liste_existe(pos_brin_1, 3)
pos_seq2b = extraire_si_pos_dans_liste_existe(pos_brin_2, 2) + ", " + extraire_si_pos_dans_liste_existe(pos_brin_2, 3)
type_hairpin_b = extraire_si_pos_dans_liste_existe(type_hpin, 1)
type_bturn_b = extraire_si_pos_dans_liste_existe(type_b_turn, 1)
zone_rama_b = extraire_si_pos_dans_liste_existe(rama, 1)

print(f"Les positions de début et de fin du brin 1 sont {pos_seq1b}")
print(f"La séquence du brin 1 est {seq1b}")
print(f"Les positions de début et de fin du brin 2 sont {pos_seq2b}")
print(f"La séquence du brin 2 est {seq2b}")
print(f"La séquence du b-hairpin est {hairpin_b}")
print(f"Le type de b-hairpin est {type_hairpin_b}")
print(f"Le type de b-turn est {type_bturn_b}")
print(f"La zone de Ramachandrant est {zone_rama_b}")
print(  )
print( )
print( )

#Vérifier si les séquences du fichier hpin existent dans la protéine
seq_proteine = ''.join(residu)
#pour seq brin 1
if seq1b in seq_proteine:
    print(f"Le brin {seq1b} existe dans la protéine")
else:
    print(f"Le brin {seq1b} n'est pas dans la protéine")

#pour seq brin 2
if seq2b in seq_proteine:
    print(f"Le brin {seq2b} existe dans la protéine")
else:
    print(f"Le brin {seq2b} n'est pas dans la protéine")

#pour seq hairpin
if hairpin_b in seq_proteine:
    print(f"Le brin {hairpin_b} existe dans la protéine")
else:
    print(f"Le brin {hairpin_b} n'est pas dans la protéine")

print()

#Pour le troisième feuillet = bloc c
seq1c = extraire_si_pos_dans_liste_existe(seq_brin_1, 2)
seq1c = seq1c.replace(" ","")
hairpin_c = extraire_si_pos_dans_liste_existe(hairpin, 2)
seq2c = extraire_si_pos_dans_liste_existe(seq_brin_2, 2)
pos_seq1c = extraire_si_pos_dans_liste_existe(pos_brin_1, 4) + ", " + extraire_si_pos_dans_liste_existe(pos_brin_1, 5)
pos_seq2c = extraire_si_pos_dans_liste_existe(pos_brin_2, 4) + ", " + extraire_si_pos_dans_liste_existe(pos_brin_2, 5)
type_hairpin_c = extraire_si_pos_dans_liste_existe(type_hpin, 2)
type_bturn_c = extraire_si_pos_dans_liste_existe(type_b_turn, 2)
zone_rama_c = extraire_si_pos_dans_liste_existe(rama, 2)

print(f"Les positions de début et de fin du brin 1 sont {pos_seq1c}")
print(f"La séquence du brin 1 est {seq1c}")
print(f"Les positions de début et de fin du brin 2 sont {pos_seq2c}")
print(f"La séquence du brin 2 est {seq2c}")
print(f"La séquence du b-hairpin est {hairpin_c}")
print(f"Le type de b-hairpin est {type_hairpin_c}")
print(f"Le type de b-turn est {type_bturn_c}")
print(f"La zone de Ramachandrant est {zone_rama_c}")
print(  )
print( )
print( )

#Vérifier si les séquences du fichier hpin existent dans la protéine
seq_proteine = ''.join(residu)
#pour seq brin 1
if seq1c in seq_proteine:
    print(f"Le brin {seq1c} existe dans la protéine")
else:
    print(f"Le brin {seq1c} n'est pas dans la protéine")

#pour seq brin 2
if seq2c in seq_proteine:
    print(f"Le brin {seq2c} existe dans la protéine")
else:
    print(f"Le brin {seq2c} n'est pas dans la protéine")

#pour seq hairpin
if hairpin_c in seq_proteine:
    print(f"Le brin {hairpin_c} existe dans la protéine")
else:
    print(f"Le brin {hairpin_c} n'est pas dans la protéine")

print()