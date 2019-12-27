# Open Chromatin and Methylation

### Лев Мазаев | мАДБМ18

---

### Этап 1: поиск данных на ENCODE

#### [CTCF ChIP-Seq](https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=TF+ChIP-seq&assembly=GRCh38&target.label=CTCF&biosample_ontology.term_name=A549) 

![scr1](/home/leo/ME/HSE/NGS/7/scr1.png)

Взял эксперимент [**ENCSR876UZD**](https://www.encodeproject.org/experiments/ENCSR876UZD/) - *Homo sapiens* A549 treated with 100 nM dexamethasone for 5 hours. Файлы bed narrowPeak: **ENCFF702GWV** и **ENCFF921XYL**.

#### [H3K4me3 ChIP-Seq](https://www.encodeproject.org/search/?type=Experiment&status=released&assembly=GRCh38&biosample_ontology.term_name=A549&target.label=H3K4me3&assay_title=Histone+ChIP-seq) 

![scr2](/home/leo/ME/HSE/NGS/7/scr2.png)

В задании указано, что нужно взять broadPeak для ChIP-Seq, но не во всех экспериментах есть такие данные. 
Подходящий эксперимент -  [**ENCSR000DPD**](https://www.encodeproject.org/experiments/ENCSR000DPD/) - *Homo sapiens* A549. Файлы bed broadPeak: **ENCFF001WUL** и **ENCFF001WUM**. Как выяснилось позднее, сами данные размечены по геному hg19, то есть нужно будет перевести их в hg38, используя [hgLiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver). 

#### [ATAC-Seq](https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq&assembly=GRCh38&biosample_ontology.term_name=A549) 

![scr3](/home/leo/ME/HSE/NGS/7/scr3.png)

Взял эксперимент [**ENCSR074AHH**](https://www.encodeproject.org/experiments/ENCSR074AHH/) - Homo sapiens A549 treated with 100 nM dexamethasone for 8 hours. Файлы bed narrowPeak: **ENCFF846JAO**, **ENCFF153HNH** и **ENCFF198GUB**.

#### [BS-Seq](https://www.encodeproject.org/search/?type=Experiment&assembly=GRCh38&assay_title=WGBS&biosample_ontology.term_name=A549)

![scr4](/home/leo/ME/HSE/NGS/7/scr4.png)

Тут только один эксперимент: [**ENCSR481JIW**](https://www.encodeproject.org/experiments/ENCSR481JIW/) - whole-genome shotgun bisulfite sequencing (WGBS) - *Homo sapiens* A549. Файлы methylation state at CpG: **ENCFF005TID** и **ENCFF003JVR**.

Скачивание производилось с помощью следующей команды:

```bash
xargs -L 1 curl -O -L < files.txt
```

### Этап 2: bed -> bigWig

Далее, распакуем все файлы с помощью `gunzip`, добавим в названия тип эксперимента и переведём их в формат bedGraph. То есть, должно остаться только 4 колонки: chr, start, end, score. Также попутно удалим все неосновные хромосомы, используя регулярное выражение  `^chr[0-9YX][0-9]?$` для первой колонки:

```bash
for file in `ls -1 ./bed`
do
        awk '$1 ~ /^chr[0-9YX][0-9]?$/ {print $1"\t"$2"\t"$3"\t"$5}' ./bed/$file \
                | sort -k1,1 -k2,2n > ./bedGraph/${file::-3}bedGraph
done
```

bedGraph для H3K4me3 реплик переразметим на геном hg38 используя [hgLiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver) и заново отсортируем (ещё пришлось объединить по 5-10 пересекающихся регионов в каждом файле при помощи bedtools merge). Теперь переведём все bedGraph в формат bigWig:

```bash
fetchChromSizes hg38 > hg38.chrom.sizes # и очистим от лишних хромосом

for file in `ls -1 ./bedGraph`
do
	./bedGraphToBigWig ./bedGraph/$file hg38.chrom.sizes ./bigWig/${file::-8}bigWig
done
```

Результат:

![scr5](/home/leo/ME/HSE/NGS/7/scr5.png)

### Этап 3: Графики обогащения

Уберём bigWig для ATAC-Seq из соответствующей папки и построим графики обогащения для всех данных  вокруг пиков первой реплики ATAC-Seq:

```bash
computeMatrix reference-point -S ./bigWig/* -R ./bed/atac_ENCFF153HNH.bed \
--referencePoint center -a 2000 -b 2000 -out me153.tab.gz -p 8
plotHeatmap -m me153.tab.gz -out me153.png --heatmapHeight 15 --heatmapWidth 7 --colorMap jet --sortRegion ascend --regionsLabel 'ATAC-Seq r153'
```

Вот что получилось:

![fail](/home/leo/ME/HSE/NGS/7/fail.jpg)

Как видим, получилось ничего. В результате гугления, выяснилось, что когда геном очень неплотно покрыт сигналом, как в нашем случае (кроме данных WGBS), нужно использовать дополнительную опцию `--missingDataAsZero`. Во всех файлах, кроме WGBS очень мало записей:

![scr6](/home/leo/ME/HSE/NGS/7/scr6.png)

Повторный запуск с новой опцией:

```bash
computeMatrix reference-point -S ./bigWig/* -R ./bed/atac_ENCFF153HNH.bed --referencePoint center -a 2000 -b 2000 -out me153.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m me153.tab.gz -out me153.png --heatmapHeight 15 --heatmapWidth 7 --regionsLabel 'ATAC-Seq r153 regions' --colorMap jet --sortRegion ascend
```

![first](/home/leo/ME/HSE/NGS/7/first.png)

Вот теперь, какой-никакой результат есть. Благодаря опции `--sortRegions ascend` регионы сортированы **по возрастанию среднего значения**. Первые две колонки - Bisulfite Sequencing для CpG метилирования, наблюдается слабый сигнал для всех пиков ATAC-Seq (целиком по каждому пику), и он везде примерно одинаковый. Видимо у этих данных слабое покрытие. Вторые две колонки - Chip-Seq транскрипционного фактора CTCF, примерно для трети ATAC-Seq регионов сигнал очень существенный в центре, но при этом узкий. При этом в тех регионах, где сигнал начинает обрываться, видна и некоторая раздельная полоса в данных BS. Видимо это из-за сортировки регионов: после раздела BS начинает влиять на среднее по региону, чем ChIP-Seq для CTCF. Для данных H3K4me3 - не видно ничего, так как score всех регионов в этих broadPeak-файлах - нулевой. 

![scr7](/home/leo/ME/HSE/NGS/7/scr7.png)

Поэтому эти данные просто не подходят для такого анализа, попробуем их заменить на что-нибудь другое, правда придется отказаться от broadPeak (таких данных очень мало и ENCODE обозначает их как плохие по качеству, ещё они все размечены по hg19) в пользу narrowPeak.

Скачаем [**ENCSR524UOX**](https://www.encodeproject.org/experiments/ENCSR524UOX/) - *Homo sapiens* A549 treated with 100 nM dexamethasone for 5 hours, H3K4me3 ChIP-seq. Файлы bed narrowPeak: **ENCFF279EVH**, **ENCFF421IEC** и **ENCFF789KGD**. Итого, получилось следующее:

![scr8](/home/leo/ME/HSE/NGS/7/scr8.png)

Будем запускать все эти файлы против каждой реплики ATAC-Seq:

```bash
computeMatrix reference-point -S ./bigWig/* -R ./atac/atac_ENCFF153HNH.bed --referencePoint center -a 2000 -b 2000 -out me153.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m me153.tab.gz -out me153.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'ATAC-Seq r153 regions' --colorMap jet --sortRegion ascend
```

![me153](/home/leo/ME/HSE/NGS/7/me153.png)

Попробуем сделать картинку без сортировки регионов по значению:

```bash
plotHeatmap -m me153.tab.gz -out me153_unsorted.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'ATAC-Seq r153 regions' --colorMap jet --sortRegion no
```

![me153_unsorted](/home/leo/ME/HSE/NGS/7/me153_unsorted.png)

Теперь у H3K4me3 настолько сильный сигнал, что вероятно из-за него не видно BS. Также этот сигнал довольной широкий, простирается почти целиком на все 4000-нуклеотидное окно с центром в пике ATAC-Seq, в то время как сигнал CTCF довольно узкий и сосредоточен только в области пика. Ещё CTCF и H3K4me3 не скоррелированы по регионам - из-за их сортировки по среднему. Тут может быть и биологическая причина: когда метилирование слишком сильное, CTCF не так часто присоединяется к ДНК, а при умеренном метилировании - максимально часто. Другая реплика:

```bash
computeMatrix reference-point -S ./bigWig/* -R ./atac/atac_ENCFF198GUB.bed --referencePoint center -a 2000 -b 2000 -out me198.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m me198.tab.gz -out me198.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'ATAC-Seq r198 regions' --colorMap jet --sortRegion ascend
```

![me198](/home/leo/ME/HSE/NGS/7/me198.png)

Такой же результат. Последняя ATAC-Seq реплика:

```bash
computeMatrix reference-point -S ./bigWig/* -R ./atac/atac_ENCFF846JAO.bed --referencePoint center -a 2000 -b 2000 -out me846.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m me846.tab.gz -out me846.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'ATAC-Seq r846 regions' --colorMap jet --sortRegion ascend
```

![me846](/home/leo/ME/HSE/NGS/7/me846.png)

От реплики к реплике результат остается один и тот же. Что в целом можно сказать:

- У Bisulfite Sequencing видимо очень слабое покрытие, сигнал был виден только до того, как мы обновили H3K4me3 данные. Сам сигнал равномерно покрывает всё 4000-нуклеотидное окно с центром в пике ATAC-Seq. Т. е. с биологической точки зрения наблюдается CpG-метилирование всего окна вокруг пика ATAC-Seq.

- CTCF данные дают узкий пик по центру окна. Так как регионы у нас сортированы по среднему значению (и видимо по всем входным трекам), то, глядя на картинки, можно сделать вывод, что для наиболее множественного соединения CTCF с  ДНК излишнее метилирование H3K4 вредно. Узкий сигнал CTCF расположен там же, где и пик ATAC-Seq, что логично, так как CTCF садится на ДНК именно в середине участка открытого хроматина, т. е. на пике ATAC-Seq. 

- Первые данные по метилированию гистоновых концов оказались очень низкого качества, score у всех регионов был нулевой, поэтому пришлось взять другие. Что видим на этих данных: в хороших по среднему значению регионах сигнал очень высокий по всему окну целиком, с уменьшением среднего значения ширина сигнала начинает падать, амплитуда тоже. Биологически это означает, что есть регионы открытого хроматина, где почти все гистоновые концы метилированы, но есть и такие регионы, где модификации гистоновых концов наблюдаются только в центре, и этого достаточно для доступа к ДНК. При таком "умеренном" сценарии также увеличивается активность CTCF . 

Посмотрим, что получится если вычислить матрицу по всем репликам:

```bash
computeMatrix scale-regions -S ./bigWig/* -R ./atac/* -b 1000 -a 1000 -out all_scaled.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m all_scaled.tab.gz -out all_scaled.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'ATAC-Seq scaled' --colorMap jet --sortRegion ascend
```

К сожалению не удалось:

> MemoryError: Unable to allocate array with shape (221515, 2100) and data type float64

Попробуем ещё 3 bigWig из данных CTCF (**ENCFF027FJN**, **ENCFF253ADS**, **ENCFF441JCU** - это именно bigWig, скачанные с [**ENCSR876UZD**](https://www.encodeproject.org/experiments/ENCSR876UZD/), далее их использовать не будем) против первой реплики:

```bash
computeMatrix reference-point -S ./CTCFdownloadedbigWigs/* -R ./atac/atac_ENCFF153HNH.bed --referencePoint center -a 2000 -b 2000 -out exp_bw.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m exp_bw.tab.gz -out exp_bw.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'ATAC-Seq r153 regions' --colorMap jet --sortRegion ascend
```

![exp_bw](/home/leo/ME/HSE/NGS/7/exp_bw.png)

Здесь результат вполне согласуется с тем, что мы видели до того: узкие линии, совпадающие по положению с пиком ATAC-Seq. При это наблюдается слабый сигнал и вокруг основного: видимо причина в том, что эти bigWig сделаны напрямую с BAM-ов, то есть это просто покрытие ридами, в то время как наши bigWig-и сделаны из данных после процедуры peak calling.

### Этап 4: Профили всех сигналов вокруг TSS человека

Воспользуемся сервисом Biomart, чтобы скачать позиции TSS:

![scr9](/home/leo/ME/HSE/NGS/7/scr9.png)

Таким образом, есть bed-файл без стренда, соответственно первое число всегда будет браться в качестве начала региона. Дополнение названий хромосом для соответствия аннотации UCSC, а также `bedtools merge`, чтобы сократить число регионов (те, что пересекаются - скорее всего один и тот же ген с разными TSS):

```bash
awk '{a="chr"$0; print a}' mart_export.txt | sort -k1,1 -k2,2n > ~/BioData/NGS_HW7/tss_raw.bed
bedtools merge tss_raw.bed > tss.bed
```

Также Теперь к уже использованным ранее bigWig, которые мы получили из разных bed-файлов, добавим ещё 3 новых с ATAC-Seq данных:

![scr10](/home/leo/ME/HSE/NGS/7/scr10.png)

И сделаем матрицу против TSS:

```bash
computeMatrix reference-point -S ./bigWig/* -R tss.bed --referencePoint TSS -a 2000 -b 2000 -out tss.tab.gz -p 8 --missingDataAsZero
plotHeatmap -m tss.tab.gz -out tss.png --heatmapHeight 20 --heatmapWidth 7 --regionsLabel 'TSS' --colorMap jet --sortRegion ascend
```

![tss](/home/leo/ME/HSE/NGS/7/tss.png)

В данном случае центр каждого окна - точка начала транскрипции. Все данные - ATAC-Seq, CTCF ChIP-Seq и H3K4me3 ChIP-Seq центрируются на этой же точке. Это говори о том, что это открытый хроматин, что здесь CTCF, и что присутствует метилирование концов гистонов. Данным BS видимо опять не хватает интенсивности.