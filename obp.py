#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:20:40 2019

@author: Alex
"""

import random


# -------------------------------------------------------------------------------------------------------------#

# FONCTION

def alphabet(type):
    """
        args : La fonction prend en argument le type d'alphabet
        return: La fonction retourne le type de lettre qui compose la séquence d'ADN
    """
    if type == "nucleic":
        return "ATGC"


def randseq(num, alpha):
    """
        args : La fonction prend en argument la taille de la séquence et l'alphabet utilisé
        return: La fonction retourne une sequence d'ADN aléatoire de la taille choisie
    """
    seq = ""
    for alphacount in range(num):
        irand = random.randint(0, len(alpha) - 1)
        seq = seq + alpha[irand]
    return seq

def cut(seq):
    """
        args : La fonction prend en argument la séquence d'ADN
        return: La fonction retourne l'ensemble des reads d'une taille choisie par l'utilisateur
                avec une liste de qualité associée
    """
    while True:
        bError = False
        choixcut = input("\nSelectionner une taille de read entre 20 et 100 :\n")
        try:
            choixcut = int(choixcut)
        except ValueError:
            print("ERREUR : Valeur incorrect")
            bError = True

        if not bError:
            if 20 <= choixcut <= 100:
                break
            else:
                print("ERREUR : La valeur n'est pas comprise entre 20 et 100")
    listcutseq = []
    listqualGlob = []
    for icountseq in range(0, len(seq), choixcut):
        listcutseq.append(seq[icountseq: icountseq + choixcut])
        listqual = []
        for icountqual in range(choixcut):
            icountqual = round(random.uniform(0.1, 1), 2)
            listqual.append(icountqual)
        listqualGlob.append(listqual)

    for i in range(len(listcutseq)):
        print(listcutseq[i], listqualGlob[i])
    return listcutseq, listqualGlob


def liste_adaptateur():
    """
        args : La fonction prend aucun argument
        return: La fonction retourne une liste d'adaptateur choisi par l'utilisateur
    """
    listeadaptateur = []
    nb_adaptateur = 0
    while True:
        nb_adaptateur = input("\nCombien d'adaptateur ?\n")
        try:
            nb_adaptateur = int(nb_adaptateur)
            break
        except ValueError:
            print("ERREUR : Valeur incorrect")

    print("Taper vos adapteurs :")
    for icountadap in range(nb_adaptateur):
        adaptateur = input()
        adaptateur = adaptateur.upper()
        listeadaptateur.append(adaptateur)
    return listeadaptateur


def nombre_adaptateur(listcutseq, adaptateur):
    """
        args : La fonction prend en argument la liste des reads et la liste des adaptateurs
        return: La fonction retourne le nombre d'adaptateur present dans le readset
    """
    count_adaptateur = 0
    for icountlistseq in range(len(listcutseq)):
        for icountadapt in range(len(adaptateur)):
            if adaptateur[icountadapt] in listcutseq[icountlistseq]:
                count_adaptateur = count_adaptateur + 1
    print("\nLe nombre total d'adaptateur est :", count_adaptateur)
    return count_adaptateur


def filtre_adaptateurs(listcutseq, listqual, adaptateur):
    """
        args : La fonction prend en argument la liste des reads, la liste des qualités associées aux reads et
               la liste des adaptateurs
        return: La fonction retourne la liste des reads et des qualités après la supression des adaptateurs présents
    """
    for icountlistseq in range(len(listcutseq)):
        for icountadapt in range(len(adaptateur)):
            for i in range(2):
                if adaptateur[icountadapt] in listcutseq[icountlistseq]:
                    adaptateurLen = len(adaptateur[icountadapt])
                    index = listcutseq[icountlistseq].index(adaptateur[icountadapt])
                    listcutseq[icountlistseq] = listcutseq[icountlistseq].replace(adaptateur[icountadapt], "")
                    for i in range(adaptateurLen):
                        listqual[icountlistseq].pop(index)
                adaptateur[icountadapt] = adaptateur[icountadapt][::-1]

    for i in range(len(listcutseq)):
        print(listcutseq[i], listqual[i])
    return listcutseq, listqual


def filtre_extrémité(listcutseq, listqual):
    """
        args : La fonction prend en argument la liste des reads et la liste des qualités associées aux reads
        return: La fonction retourne la liste des reads et des qualités après la supression des qualités des extrmités
                inférieur au seuil défini par l'utilisateur
    """
    while True:
        bError = False
        seuil = input("\nEntrez le seuil minimal de qualité compris entre 0 et 1 :\n")
        try:
            seuil = float(seuil)
        except ValueError:
            print("ERREUR : Valeur incorrect")
            bError = True

        if not bError:
            if 0.0 <= seuil <= 1.0:
                break
            else:
                print("ERREUR : La valeur n'est pas comprise entre 0 et 1")
    lsQual = []
    lsNu = []

    for icountlistqual in range(len(listqual)):
        lsSaveQ = []
        lsSaveN = []
        testSup = False
        for icountlistqual2 in range(len(listqual[icountlistqual])):
            if testSup:
                lsSaveQ.append(listqual[icountlistqual][icountlistqual2])
                lsSaveN.append(listcutseq[icountlistqual][icountlistqual2])
            elif listqual[icountlistqual][icountlistqual2] >= seuil:
                lsSaveQ.append(listqual[icountlistqual][icountlistqual2])
                lsSaveN.append(listcutseq[icountlistqual][icountlistqual2])
                testSup = True
        lsNu.append(lsSaveN)
        lsQual.append(lsSaveQ)

    for icountlistqual in range(len(lsQual)):
        for icountlistqual2 in range(len(lsQual[icountlistqual]) - 1, 0, -1):
            if lsQual[icountlistqual][icountlistqual2] <= seuil:
                del lsQual[icountlistqual][icountlistqual2]
                del lsNu[icountlistqual][icountlistqual2]
            else:
                break

    for i in range(len(lsNu)):
        for j in range(len(lsNu[i])):
            print(lsNu[i][j], end='')
        print(" ", lsQual[i])

    return lsNu, lsQual


def filtre_moyenne(listcutseq, listqual):
    """
        args : La fonction prend en argument la liste des reads et la liste des qualités associées aux reads
        return: La fonction retourne la liste des reads et des qualités après la supression des reads avec une qualité
                moyenne inférieur au seuil défini par l'utilisateur
    """
    while True:
        bError = False
        seuilMOY = input("\nEntrez le seuil moyen de qualité compris entre 0 et 1 :\n")
        try:
            seuilMOY = float(seuilMOY)
        except ValueError:
            print("ERREUR : Valeur incorrect")
            bError = True

        if not bError:
            if 0.0 <= seuilMOY <= 1.0:
                break
            else:
                print("ERREUR : La valeur n'est pas comprise entre 0 et 1")

    lsQual = []
    lsNu = []

    for i in range(len(listqual)):
        somme = 0
        for j in range(len(listqual[i])):
            somme = somme + listqual[i][j]

        moyenne = somme / len(listqual[i])

        if moyenne >= seuilMOY:
            lsQual.append(listqual[i])
            lsNu.append(listcutseq[i])

    for i in range(len(lsNu)):
        for j in range(len(lsNu[i])):
            print(lsNu[i][j], end='')
        print(" ", lsQual[i])

    return lsNu, lsQual


def filtre_read(listcutSeq):
    """
        args : La fonction prend en argument la liste des reads
        return: La fonction retourne le ou les reads qui produisent le chevauchement maximal
                par rapport au motif donné par l'utilisateur
    """
    readInput = input("\nEntrez le read a comparer : \n")
    readInput = readInput.upper()

    nbPointTotal = []
    lsNucMax = []
    for i in range(len(listcutSeq)):
        nbPoint = 0
        for j in range(len(readInput)):
            if readInput[j] == listcutSeq[i][j]:
                nbPoint += 1
        for k in range(len(readInput)):
            if readInput[k] == listcutSeq[i][len(listcutSeq[i]) - k - 1]:
                nbPoint += 1

        nbPointTotal.append(nbPoint)

    maxPoint = max(nbPointTotal)

    for i in range(len(listcutSeq)):
        if nbPointTotal[i] == maxPoint:
            lsNucMax.append(listcutSeq[i])

    print("Le(s) read(s) qui produit/produisent le meilleur chevauchement avec {} ressemblance(s) : \n".format(maxPoint))
    for i in range(len(lsNucMax)):
        for j in range(len(lsNucMax[i])):
            print(lsNucMax[i][j], end='')
        print("")

    return lsNucMax

def flatten(listcutSeq):
    """
        args : La fonction prend en argument la liste des reads
        return: La fonction retourne la séquence génomique après les étapes de filtration
    """
    strFinal = ""
    for i in range(len(listcutSeq)):
        for j in range(len(listcutSeq[i])):
            strFinal = strFinal + listcutSeq[i][j]

    print("\nLa séquence génomique généré avec le readset filtré est :")
    print(strFinal)
    print("Elle a une taille de {}".format(len(strFinal)))

    return strFinal


# -------------------------------------------------------------------------------------------------------------#

# MAIN

render = False  # Mettre a True pour voir les aides des fonctions

alphabet = alphabet("nucleic")

seq = randseq(1000, alphabet)
if render:
    help(randseq)
for i in range(len(seq)):
    print(seq[i], end='')

listcutSeq, listqual = cut(seq)
if render:
    help(cut)

adaptateur = liste_adaptateur()
if render:
    help(liste_adaptateur)

nombre_adaptateur(listcutSeq, adaptateur)
if render:
    help(nombre_adaptateur)

listcutSeq, listqual = filtre_adaptateurs(listcutSeq, listqual, adaptateur)
if render:
    help(filtre_adaptateurs)

listcutSeq, listqual = filtre_extrémité(listcutSeq, listqual)
if render:
    help(filtre_extrémité)

listcutSeq, listqual = filtre_moyenne(listcutSeq, listqual)
if render:
    help(filtre_moyenne)

lsNucMax = filtre_read(listcutSeq)
if render:
    help(filtre_read)

SeqFinal = flatten(listcutSeq)
if render:
    help(flatten)
