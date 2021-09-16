#!/bin/bash
#Clearen das Qc-Data verzeichnisses um kollisionen mit alten Daten zu vermeiden
cd /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor
rm *.txt


####Daten aus Hardklor output auslesen

#lesen der 5 Spalte (Base Isotope Peak) der 1. (als indentifier welche Zeilen Scans waren) sowie der Zeilennummer (als späterer Index) der .hk Datei | und löschen der 4. Nachkommastelle
awk '{print NR, " " $1 " " $3 " " $5 " " $4}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/centroided-sample-human-orbi2-V1-S3.hk  > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt

#Erstellen eines Arrays, dass alle Zeilenzahlen der Scanzeilen enthällt
declare -a scanlines
allscanlines=$(awk '/S/ {print NR}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt)

scanlines=(${allscanlines// / })

## Die RT_scans Datei wurde aus dem Output von Decon2LS erstellt, da dort die RT in Sekunden angegeben ist und ich so Rundungsfehlern aus dem weg gehe. Allerdings muss dort noch die Zeile 1500
## ans Ende angefügt werden, da Hardklor einen Scan mehr gefunden hat als Decon2LS (501 anstatt nur 500)

#Erstellung eines Files mit allen Scanlines

printf "%s\n" "${scanlines[@]}" > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/scaninfo.txt

# Mergen der Scanlines mit den dazugehörigen  RT Werte -> File mit Lines die Den zeilen Index jedes Scans mit der dazugehörigen RT des Scans enthällt
paste -d ' ' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/scaninfo.txt /buffer/ag_bsc/pmsb_2021/Rieckert/data/RT_scans.txt > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/scans_index_RT.txt


#Printen der Base Isotope Peak werte  HK Zeilen zusammen mit ihrem Zeilen Index, löschen der S Zeilen die aus den Scan-Zeilen stammen. Diese stellen keine Peaks dar
awk '{print $1 " " $2 " " $4 " " $5}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt | sed '/S/d' > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHK.txt
awk '{print $1 " " $3 " " $4}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHK.txt >  /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKFinal.txt

#Löschen der zwischenspeicherdatein
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt
#rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHK.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/scaninfo.txt

echo "Array Erstellung"
declare -A peakArray
declare -A peakArray2
declare -A scanArray
declare -A peakRTArray
declare -a peakIndex
declare -a scanIndex
declare -A minPeakArray
declare -A maxPeakArray
trenner="/"



#Aufbau des peakArrays zum Sammeln aller Zeilenindexe der Peaks zusammen mit den entsprechenden MZ Values Es handelt sich also um ein Dictonary (Array mit Keys=Zeilenindex und Values=MZ values)
while IFS= read -r line
do
        tmpKey=$(echo $line | cut -d' ' -f1)
        tmpVal=$(echo $line | cut -d' ' -f2)
        peakArray["$tmpKey"]="$tmpVal"
        peakIndex+=($tmpKey)
#       echo "${peakArray[$tmpKey]}"

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKFinal.txt
#echo "${peakArray[2]}"
#echo "${peakIndex[1]}"


while IFS= read -r line
do
        tmpKey=$(echo $line | cut -d' ' -f1)
        tmpVal=$(echo $line | cut -d' ' -f3)
        peakArray2["$tmpKey"]="$tmpVal"
        #peakIndex+=($tmpKey)
#       echo "${peakArray2[$tmpKey]}"

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKFinal.txt



#Aufbau des scanArrays aller Zeilenindexe der Scans zusammen mit den entsprechenden RT Values Es handelt sich also um ein Dictonary (Array mit Keys=Zeilenindex und Values=RT values)
while IFS= read -r line
do
        tmpKey=$(echo $line | cut -d' ' -f1)
        tmpVal=$(echo $line | cut -d' ' -f2)
        scanArray["$tmpKey"]="$tmpVal"
        scanIndex+=($tmpKey)

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/scans_index_RT.txt
#echo "${scanArray[1]}"






####Erstellung der finalen Peakdatei

#Index Variablen für bessere lesbarkeit der Schleifenköpfe
numberScans=${#scanIndex[@]}
numberScansKleiner=$((numberScans-1))

#echo "Erstellung der finalen Peakdatei."

#Loops zum erstellen der Finalen Peaks Datei die zu jedem Peak auch die RT mit aufführt
for i in "${peakIndex[@]}"; do
        #merken des Zeilenindex des im moment betrachteten Peaks
        PeakInd=$i
        #Nun wird für den Peak gecheckt in welchem Scan er drin liegt
        for ((index=0; index <  $numberScans; index++)); do
                #Merken des Zeilen index des kleineren Scans
                ScanInd=${scanIndex[index]}
                #Falls wir im Letzten Scan angekommen sind wird überprüft ob der zeilenindex des Peaks größer ist als der des Scans, falls ja werden die Peaks mit der RT des letzten scans ausgegeben
                if [ $((index)) -eq $((numberScansKleiner)) ] && [ $((PeakInd)) -gt $((ScanInd)) ];then
                         echo "${peakArray[$PeakInd]} ${scanArray[$ScanInd]}" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filtPeaksRT.txt
                else
                        #Falls nicht im letzten Scan wird überprüft ob der Peakindex zwischen dem des betrachteten und des nächsten Scans liegt, falls ja wird der Peak mit der RT des betrachteten Scans
                        #ausgegeben
                        nextScanInd=${scanIndex[index+1]}
                        if [ $((PeakInd)) -gt $((ScanInd)) ] && [ $((PeakInd)) -lt $((nextScanInd)) ];then
                                echo "${peakArray[$PeakInd]} ${scanArray[$ScanInd]} ${peakArray2[$PeakInd]}" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filtPeaksRT.txt
                        fi
                fi
        done
done

