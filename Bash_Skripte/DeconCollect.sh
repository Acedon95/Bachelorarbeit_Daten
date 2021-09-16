#!/bin/bash

awk -F "," '{print $1 " " $4}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/raw-sample-human-orbi2_isos.csv  > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/index_peakTMP.txt
tail -n +2 "/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/index_peakTMP.txt" > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/index_peak.txt
#Lesen der Scan Info |

awk -F "," '{print $1 " " $2}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/raw-sample-human-orbi2_scans-Decon2LS.csv > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/RT_scansTMP.txt
tail -n +2 "/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/RT_scansTMP.txt" > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/RT_scans.txt

awk -F '\t' '{print $2 " " $3 " " $4}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/raw-sample-human-orbi2_peaks-Decon2LS.txt > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/Int_scansTMP.txt
tail -n +2 "/buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/Int_scansTMP.txt" > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/Int_scans.txt

echo "Array Erstellung"
declare -A peakArray
declare -A scanArray
declare -A intArray
declare -a peakIndex
declare -a scanIndex


trenner="/"
tmppeakCout=0
#Aufbau des peakArrays zum Sammeln aller Scanindexe der Peaks zusammen mit den entsprechenden MZ Values Es handelt sich also um ein Dictonary (Array mit Keys=Scanindex und Values=MZ values)
while IFS= read -r line
do
        tmppeakCout=$((tmppeakCout+1))
        tmpKey=$(echo $line | cut -d' ' -f1)
        tmpVal=$(echo $line | cut -d' ' -f2)
        peakArray["$tmpKey$trenner$tmppeakCout"]="$tmpVal"
        peakIndex+=("$tmpKey$trenner$tmppeakCout")

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/index_peak.txt

#Aufbau des scanArrays aller Zeilenindexe der Scans zusammen mit den entsprechenden RT Values Es handelt sich also um ein Dictonary (Array mit Keys=Zeilenindex und Values=RT values)
while IFS= read -r line
do
        tmpKey=$(echo $line | cut -d' ' -f1)
        tmpVal=$(echo $line | cut -d' ' -f2)
        scanArray["$tmpKey"]="$tmpVal"
        scanIndex+=($tmpKey)

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/RT_scans.txt

while IFS= read -r line
do
        tmpKey1=$(echo $line | cut -d' ' -f1)
	tmpKey2=$(echo $line | cut -d' ' -f2)
        tmpVal=$(echo $line | cut -d' ' -f3)
        intArray["$tmpKey1$trenner$tmpKey2"]="$tmpVal"
	
	#echo "SCannum: $tmpKey1 MZ: $tmpKey2 Int: $tmpVal Key: $tmpKey1$trenner$tmpKey2 ArrayWert: ${intArray[$tmpKey1$trenner$tmpKey2]}"

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/Int_scans.txt



numberScans=${#scanIndex[@]}

echo "Erstellung der finalen Peakdatei."


#Loops zum erstellen der Finalen Peaks Datei die zu jedem Peak auch die RT mit aufführt
for i in "${peakIndex[@]}"; do
        #merken des Zeilenindex des im moment betrachteten Peaks
        PeakInd=$i
        PeakIndMod=${PeakInd%/*}
	PeakMz=${peakArray[$PeakInd]}
        #Nun wird für den Peak gecheckt in welchem Scan er drin liegt
        for ((index=0; index <  $numberScans; index++)); do
                #Merken des Zeilen index des kleineren Scans
                ScanInd=${scanIndex[index]}
                #Falls wir im Letzten Scan angekommen sind wird überprüft ob der zeilenindex des Peaks größer ist als der des Scans, falls ja werden die Peaks mit der RT des letzten scans ausgegeben
                if [ $((PeakIndMod)) -eq $((ScanInd)) ];then
                        #Falls nicht im letzten Scan wird überprüft ob der Peakindex  dem des betrachteten Scans gleicht, falls ja wird der Peak mit der RT des betrachteten Scans
                        #ausgegeben
                        echo "${scanArray[$ScanInd]} ${peakArray[$PeakInd]} ${intArray[$ScanInd$trenner$PeakMz]}" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcDecon/Decon2LS.dta2d
                fi
        done
done
