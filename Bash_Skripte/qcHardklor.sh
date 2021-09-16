#!/bin/bash
#Clearen das Qc-Data verzeichnisses um kollisionen mit alten Daten zu vermeiden
cd /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor
rm *.txt


####Daten aus Hardklor output auslesen

#lesen der 5 Spalte (Base Isotope Peak) der 1. (als indentifier welche Zeilen Scans waren) sowie der Zeilennummer (als späterer Index) der .hk Datei | und löschen der 4. Nachkommastelle
awk '{print NR, " " $1 " " $3 " " $5}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/centroided-sample-human-orbi2-V1-2.hk  > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt

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
awk '{print $1 " " $2 " " $4}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt | sed '/S/d' > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHK.txt
awk '{print $1 " " $3}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHK.txt >  /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKFinal.txt

#Löschen der zwischenspeicherdatein
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKTMP.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHK.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/scaninfo.txt



#### Auslesen der Daten aus der FeatureXML 

#Lesen der MZ Values der gefundenen Features | cut befehle zum Löschen der XML Notationen
sed -n "/MZ=/p" /buffer/ag_bsc/pmsb_2021/Rieckert/data/featureMap-sample-human-orbi2.featureXML| cut -d'=' -f6 | cut -d' ' -f1 >  /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredFeatureMZ.txt #| cut -d'=' -f6 | cut -d' ' -f1

#Die folgenden 3 zeilen dienen zum entfernen der verbliebenen " und dem löschen der Nachkommastellen ab stelle 4 um das Problem des durch Hardklor vorgenommenen Runden zu umgehen 
#(Ausbaufähig durch rundungsfkt)
sed -n "/./p" /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredFeatureMZ.txt | cut -d'.' -f1 | cut -c2- > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/vorKomma.txt
sed -n "/./p" /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredFeatureMZ.txt | cut -d'.' -f2 | sed 's/.$//' > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/nachKomma.txt
paste -d "." /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/vorKomma.txt /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/nachKomma.txt > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredFeatureMZFinal.txt

#Lesen der Convexen Hüllen nr 0 (Monoisotopische) Für jedes feature | Es werden die zwei nachfolgenden zeilen nach dem match von "convexhull nr=\"0\"" ausgegeben und Teile der XML Notation gelöscht 
# genau: <convexhull  und <pt  sowie  y="...." />
grep -ws "convexhull nr=\"0\"" /buffer/ag_bsc/pmsb_2021/Rieckert/data/featureMap-sample-human-orbi2.featureXML -A 2 |awk '{print $2}' > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0_TMP1.txt

#Entfernen des nr=
awk '{print $1}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0_TMP1.txt | sed "/nr=/d" > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0_TMP2.txt

#Entfernen des X=" und des abschließenden "
cut -d'"' -f2 /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0_TMP2.txt > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0_TMP3.txt

#Zusammenfügern der Grenzen der convexen Hüllen in eine Zeile | Erst leerzeilen entfernung mit sed und dann mergen von Zeilen 1 und 2, 3 und 4 usw mittels paste
sed '/^[[:space:]]*$/d' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0_TMP3.txt | paste -d " " - - > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0Final.txt

#Lesen der elution_profile_intensities Informationen aus dem XML sowie entfernen der XML Notation  ---> Verbleibend [...........] Intensities Liste 
sed -n "/elution_profile_intensities/p" /buffer/ag_bsc/pmsb_2021/Rieckert/data/featureMap-sample-human-orbi2.featureXML | cut -d'=' -f4 | cut -d'"' -f2  > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/elutionProfile_TMP1.txt

#Zählen der Elemente pro elution_profile_intensities Liste
awk '{print NF}' /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/elutionProfile_TMP1.txt > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/elution_profile_count.txt

#verbinden aller temporären Feature Daten zu einer großen Datei
paste -d " " /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredFeatureMZFinal.txt /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/convHull0Final.txt /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/elution_profile_count.txt > /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/finalFeatures.txt



####Einfügen der RT Werte in die Peaklisten Datei, sodass jeder Peak die RT des Scnas erhällt aus dem er stammt
echo "Array Erstellung"
declare -A peakArray
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
#	echo "${peakArray[$tmpKey]}"

done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filteredHKFinal.txt
#echo "${peakArray[2]}"
#echo "${peakIndex[1]}"

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
				echo "${peakArray[$PeakInd]} ${scanArray[$ScanInd]}" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filtPeaksRT.txt
			fi
		fi
	done
done

##Aufbau von Peakarray mit key rt+trenner (zur uniq) value mz

tmppeakCout=0
while IFS= read -r line
do
        tmppeakCout=$((tmppeakCout+1))
        tmpKey=$(echo $line | cut -d' ' -f2)
        tmpVal=$(echo $line | cut -d' ' -f1)
        #Runden auf 2 nachkommastellen nach dem ersetzen von . zu ,
        tmpVal=${tmpVal/./,}
        tmpVal=`printf "%.2f\n" "$tmpVal"`
        tmpVal=${tmpVal/,/.}
        peakRTArray["$tmpKey$trenner$tmppeakCout"]="$tmpVal"


done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/filtPeaksRT.txt

####beginn der Datenauswertung
echo "Beginn der Datenauswertung"
#Anzahl der Peaks die gefunden wurden
foundCount=0
existCount=0
abweichungMz=0.01
negAbweichCount=0
posAbweichCount=0
tmpChar="sek"
declare -A foundPeaksFeat

while IFS= read -r line
do
	echo "Start neues feature"
	#Auslesen der daten des betrachteten Features
	mzFeat=$(echo $line | cut -d' ' -f1)
	#Runden auf 2 nachkommastellen nach dem ersetzen von . zu , da printf für das runden ein , braucht
        mzFeat=${mzFeat/./,}
        mzFeat=`printf "%.2f\n" "$mzFeat"`
        #rückwandeln zu . da , einen stdin fehler in späteren abfragen wirft
        mzFeat=${mzFeat/,/.}
	rtMin=$(echo $line | cut -d' ' -f2)
	#rtMin=$((rtMin+)) HIER MUSS DIE BERECHNUNG FÜR DIE ERLAUBTE ABWEICHUNG HERREIN!!
	rtMax=$(echo $line | cut -d' ' -f3)
	peakCountFeat=$(echo $line | cut -d' ' -f4)
	#Liste der gefundenen Peaks im feature und counter für die Anzahl der gefundenen
	unset foundPeaksFeat
	foundPeaksFeatCount=0
	#zur unique haltung der Keys im Asugabe array
	controllCount=0


        for x in "${!peakRTArray[@]}";do
                controllCount=$((controllCount+1))
                #Auslesen der peak Informationen m/z und RT des Peaks
                mzPeak=${peakRTArray[${x}]}
               # mzPeak=${mzPeak/./,}
               # mzPeak=`printf "%.2f\n" "$mzPeak"`
               # mzPeak=${mzPeak/,/.}
                rtPeakTMP=$x
                rtPeak=${x%/*}
                minMzPeak=`echo "$mzPeak - $abweichungMz"|bc -l`
               maxMzPeak=`echo "$mzPeak + $abweichungMz"|bc -l`
		#Check ob Peak in der convexen Huelle 0 des features liegt  diese Variate mit bc ist nötig das es sich um Float Variablen handelt mit denen Bash nicht umgehen kann
		#bc ist eine Taschenrechner anwendung und gibt 0 oder 1 zurück wenn das statement true oder false ist
		if (( $(echo "$rtPeak > $rtMin" |bc -l) ));then
			if (( $(echo "$rtPeak < $rtMax" |bc -l) ));then
				if (( $(echo "$mzPeak == $mzFeat" |bc -l) ));then
					foundPeaksFeatCount=$((foundPeaksFeatCount+1))
					foundCount=$((foundCount+1))
					foundPeaksFeat["$rtPeak$controllCount"]="$mzPeak"
				else
	                                if (( $(echo "$minMzPeak == $mzFeat" |bc -l) ));then
	                                        foundPeaksFeatCount=$((foundPeaksFeatCount+1))
	                                        foundCount=$((foundCount+1))
						negAbweichCount=$((negAbweichCount+1))
	                                        foundPeaksFeat["$rtPeak$controllCount"]="$mzPeak"
					else
		                                if (( $(echo "$maxMzPeak == $mzFeat" |bc -l) ));then
		                                        foundPeaksFeatCount=$((foundPeaksFeatCount+1))
		                                        foundCount=$((foundCount+1))
							posAbweichCount=$((posAbweichCount+1))
		                                        foundPeaksFeat["$rtPeak$controllCount"]="$mzPeak"
		                                fi
					fi
				fi
			fi
		fi
	done

	#Erhöhen der Anzahl der Peaks die in der gesammten FeatureMap existieren
	existCount=$((existCount+peakCountFeat))
	echo "Im Feature mit m/z=$mzFeat und der convexen Huelle von $rtMin - $rtMax wurden $foundPeaksFeatCount von $peakCountFeat Peaks gefunden. Diese sind:" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/qcSummary.txt
	for x in "${!foundPeaksFeat[@]}";do
		echo "Peak mit RT= $x m/z=${foundPeaksFeat[$x]}." >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/qcSummary.txt
	done
done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/finalFeatures.txt

echo "Von $existCount in den Featuren enthaltenen Peaks wurden $foundCount erfolgreich von Hardklor gefunden. Davon $negAbweichCount mit einer M/z Abweichung von -0.01 und $posAbweichCount mit einer von +0.01" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/qcSummary.txt


falseFoundCount=0
peakCount=0
mzPeak=""
for y in "${!peakRTArray[@]}";do
        peakCount=$((peakCount+1))
        found=1
        #Auslesen der peak Informationen m/z und RT des Peaks
        rtPeakTMP=$y
        rtPeak=${y%/*}
        mzPeak=${peakRTArray[${y}]}
#        mzPeak=${mzPeak/./,}
#        mzPeak=`printf "%.2f\n" "$mzPeak"`
#        mzPeak=${mzPeak/,/.}
        minMzPeak=`echo "$mzPeak - $abweichungMz"`
        maxMzPeak=`echo "$mzPeak + $abweichungMz"`
	#Hier Fehlt der abgleich mit den erlaubten Abweichungen (evt mittels aritmetische operation innerhalb der if bedingung)
	#Check ob Peak in der convexen Huelle 0 des features liegt  diese Variate mit bc ist nötig das es sich um Float Variablen handelt mit denen Bash nicht umgehen kann
	#bc ist eine Taschenrechner anwendung und gibt 0 oder 1 zurück wenn das statement true oder false ist
	while IFS= read -r line
	do
	        #Auslesen der daten des betrachteten Features
	        mzFeat=$(echo $line | cut -d' ' -f1)
	        #Runden auf 2 nachkommastellen nach dem ersetzen von . zu , da printf für das runden ein , braucht
	        mzFeat=${mzFeat/./,}
	        mzFeat=`printf "%.2f\n" "$mzFeat"`
	        #rückwandeln zu . da , einen stdin fehler in späteren abfragen wirft
	        mzFeat=${mzFeat/,/.}
	        rtMin=$(echo $line | cut -d' ' -f2)
	        rtMax=$(echo $line | cut -d' ' -f3)

		if (( $(echo "$rtPeak > $rtMin" |bc -l) ));then
			if (( $(echo "$rtPeak < $rtMax" |bc -l) ));then
				if (( $(echo "$mzPeak == $mzFeat" |bc -l) ));then
					found=0
					falseFoundCount=$((falseFoundCount+1))
					break
				else
	                                if (( $(echo "$minMzPeak == $mzFeat" |bc -l) ));then
	                                        #foundPeaksFeatCount=$((foundPeaksFeatCount+1))
	                                        falseFoundCount=$((falseFoundCount+1))
	                                        found=0
						break
						# foundPeaksFeat["$rtPeak"]="$mzPeak"
	                                else
	                                	if (( $(echo "$maxMzPeak == $mzFeat" |bc -l) ));then
	                                	        #foundPeaksFeatCount=$((foundPeaksFeatCount+1))
	                                                falseFoundCount=$((falseFoundCount+1))
	                                	        found=0
							break
							# foundPeaksFeat["$rtPeak"]="$mzPeak"
	                                	fi
					fi
				fi
			fi
		fi

	done < /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/finalFeatures.txt
	
	if [ found != 0 ]; then
		echo "Der Peak mit m/z = $mzPeak, RT = $rtPeak konnte in keinem Feature gefunden werden" >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/falseFound.txt
	fi

done

falseFoundCount=$((peakCount-falseFoundCount))

echo "Von den $peakCount von Hardklor gefundenen Peaks waren $falseFoundCount falsch gefunden." >> /buffer/ag_bsc/pmsb_2021/Rieckert/data/qcHardklor/falseFound.txt





#Aufräumen
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/vorKomma.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/nachKomma.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/filteredFeatureMZ.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/convHull0_TMP1.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/convHull0_TMP2.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/convHull0_TMP3.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/convHull0Final.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/elutionProfile_TMP1.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/elution_profile_count.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/filteredFeatureMZFinal.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/filteredHKFinal.txt
rm /buffer/ag_bsc/pmsb_2021/Rieckert/data/qc/scans_index_RT.txt







