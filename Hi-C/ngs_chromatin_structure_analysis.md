# Chromatin Structure Analysis

### Лев Мазаев | мАДБМ18

---

### Задание 1. Yue lab browser

a:

![scr1](/home/leo/ME/HSE/NGS/8/scr1.png)

С точки зрения аннотации (что более точно) ТАД-ами ген *HBA1* находится **внутри ТАД-а**. Если приблизить изображение и поразглядывать, то

![scr2](/home/leo/ME/HSE/NGS/8/scr2.png)

*HBA1* находится близко к границе двух визуально выделяющихся треугольников, это могут быть неаннотированные ТАД-ы меньшего размера. Теперь загрузим данные тех же клеток, но с использованием VC-коррекции (b): 

![scr3](/home/leo/ME/HSE/NGS/8/scr3.png)

ТАД-ы стало видно лучше, но пропала их аннотация. Увеличим участок с геном HBA1:

![scr4](/home/leo/ME/HSE/NGS/8/scr4.png)

Видим, что он расположен на границе двух визуально выделяющихся ТАД-ов. (c) Судя по разметке DHS, рядом с этим геном очень много сайтов, чувствительных к ДНК-азе, при этом сам ген находится между такими сайтами, значит он в области закрытого хроматина, а слева и справа много открытых участков. (d) Учитывая то, что клеточная линия К562 получена из эритролейкемии (при этой болезни снижается уровень гемоглобина в крови), а HBA1 - ген альфа-глобина, то вполне ожидаемо, что активность этого гена будет снижена при этой болезни.

### Задание 2: Сравнение данных Hi-C в HiGlass

a: Первое бросающееся в глаза отличие между клеточными линиями GM12878 и K562 (шкала примерно одинаковая) - у первой затемненные участки (контакты в хроматине) значительно больше, чем у второй: визуально хорошо различается несколько крупных квадратов (ТАД-ов) слева (GM12878), которые плохо видны справа (K562), в то время как справа хорошо различимы более мелкие ТАД-ы. То есть, можно предположить, что у GM12878 несколько более плотный хроматин, обеспечивающий большее число удаленных контактов. Также на левой матрице слева вверху выделяются точки увеличенной плотности контактов в одном из ТАД-ов, они могут соответствовать петлям. У K562 такие точки в той области менее заметны. Зато у К562 выделяются петли между доменами в середине:

![scr5](/home/leo/ME/HSE/NGS/8/scr5.png)

b: Скриншот изменения параметров цветовой шкалы с целью добиться схожести:

![scr6](/home/leo/ME/HSE/NGS/8/scr6.png)

c: Не сказать, что структуры становится похожи. В целом качество экспериментов неравноценно, матрица справа кажется куда более зашумленной даже при полном ползунке. Но при этом число различимых доменов на правой картинке также выше, будто бы контрастность больше, непонятно ввиду данных ли или качества эксперимента. В целом, всё это несколько сбивает с толку при сравнении.

d, e: Открыты обе матрицы + настроено отбражение.

![scr7](/home/leo/ME/HSE/NGS/8/scr7.png)

f, g: В целом уже по картинкам выше видно, что и внутрихромосомных, так и межхромосомных контактах на матрице справа больше, особенно заметно в правом нижнем углу. Но приблизим в центре:

![scr8](/home/leo/ME/HSE/NGS/8/scr8.png)

Какая разница заметна: матрица справа теплее, значит больше контактов в целом согласно шкале. Кроме того, на ней чётче выделяются домены, относящиеся к хромосомам, обозначая большее число дальних внутрихромосомных контактов. Насчёт межхромосомных контактов - все те же петли одинаковы и слева, и справа, но общий "контактный фон" (тёплость матрицы вне хромосомных доменов) справа выше. Объяснение наблюдаемых различий следующее - слева клетки находятся в прометафазе, это одна из фаз деления клетки. В ходе неё разрушается ядерная оболочка, а хромосомы оказываются в экваториальной плоскости веретена деления. Предыдущая фаза - профаза, в ходе неё хромосомы конденсируются (становятся различимы), что очевидно приводит к уменьшению числа как дальних внутрихромосомных, так и межхромосомных контактов. Данные справа - фаза G1, период роста клетки после митоза, в это время клетка активно синтезирует мРНК и белки, при этом ДНК представляет собой плотный хроматин. 

### Задание 3: Петли, ТАДы и CTCF

a: ТАД-ы выделены синим, петли - зелёным.

![scr9](/home/leo/ME/HSE/NGS/8/scr9.png)

Наблюдаю по крайней мере 7 визуально хорошо различимых ТАД-ов, также 3 отдельных петли и одну группу петель, которая выделена овалом.

b: Сравнивая с моим выделением выше, могу сказать, что с некоторым крупными ТАД-ами было угадано, но маленькие разглядеть не представляется возможным. Группы петель также удалось выделить, но не все. 

![scr10](/home/leo/ME/HSE/NGS/8/scr10.png)

с: 

![scr11](/home/leo/ME/HSE/NGS/8/scr11.png)

Глядя на картинку, создаётся такое ощущение, что CTCF любит концентрироваться на границах ТАД-ов. Там же находятся и петли. Внутри ТАД-ов или вне их CTCF встречается более разреженно. 

### Задание 4: Настройка системы

```bas
conda create -n chr python=3.6
conda activate chr
conda install -c bioconda deeptools
conda install -c bioconda -c conda-forge cooler
conda install -c bioconda -c conda-forge hicexplorer
pip install hic2cool
```

a:

![scr12](/home/leo/ME/HSE/NGS/8/scr12.png)

b: 

![scr13](/home/leo/ME/HSE/NGS/8/scr13.png)

### Задание 5: Анализ данных Hi-C

a: Разрешение - 10 Kb, геномная сборка - dm3 (дрозофила), количество контактов - 28196414.

![scr14](/home/leo/ME/HSE/NGS/8/scr14.png)

Изменим разрешение до 20 Kb:

```bash
cooler coarsen -k 2 -o Kc167.20000.cool Kc167.10000.cool
```

Пересохраним в h5:

```bash
hicConvertFormat --matrices Kc167.20000.cool -o Kc167.20000.h5 --inputFormat cool --outputFormat h5
```

b, c: Проведём итеративную коррекцию, а именно для начала построим диагностический график:

```bash
hicCorrectMatrix diagnostic_plot -m Kc167.20000.h5 -o hic_corrected.png
```

![hic_corrected](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/hic_corrected.png)

Видим, что пороги надо выбирать следующие: -1.6 и 4, чтобы убрать выбросы. Теперь проведём коррекцию без подбора параметров и с подбором, и построим карты по порядку (исходная, наивная коррекция, коррекция с подобранными параметрами):

```bash
hicCorrectMatrix correct -m Kc167.20000.h5 --filterThreshold -10 10 -n 10 --out Kc167.corr.20000.h5
hicCorrectMatrix correct -m Kc167.20000.h5 --filterThreshold -1.6 4 -n 500 --out Kc167.corr2.20000.h5
hicPlotMatrix -m Kc167.20000.h5 -o Kc167.raw.mtx.png --log1p --clearMaskedBins --region chrX:10000000-12000000
hicPlotMatrix -m Kc167.corr.20000.h5 -o Kc167.corr.mtx.png --log1p --clearMaskedBins --region chrX:10000000-12000000
hicPlotMatrix -m Kc167.corr2.20000.h5 -o Kc167.corr2.mtx.png --log1p --clearMaskedBins --region chrX:10000000-12000000
```

![Kc167.raw.mtx](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/Kc167.raw.mtx.png)

![Kc167.corr.mtx](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/Kc167.corr.mtx.png)

![Kc167.corr2.mtx](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/Kc167.corr2.mtx.png)

Как видим, при проведении наивной коррекции часть шума с карты ушла, ТАД-ы стали чуть более заметными. Но подбор параметров коррекции ощутимо улучшил результат: шума стало сильно меньше, квадраты ТАД-ов теперь визуально прослеживаются. Собственно, теперь поищем ТАД-ы в скорректированной версии карты и посмотрим, что получилось:

```bash
hicFindTADs -m Kc167.corr2.20000.h5 --outPrefix tad1 --minDepth 60000 --maxDepth 1000000  --step 20000 --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr
hicPlotTADs --tracks tracks1.ini --region chr2L:1000000-4000000 -o tads1.png
```

![tads1](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/tads1.png)

d: ТАД-ы получились размерами от 10 Кб до порядка одной мегабазы. Параметры графика в следующем файле:

```ini
[x-axis]
where = top
fontsize = 20

[hic matrix]
file = Kc167.corr2.20000.h5
title = Hi-C
colormap = Reds
depth = 1000000
transform = log1p
x_labels = True
file_type = hic_matrix
show_masked_bins = True
scale_factor = 1

[tads]
file = tad1_domains.bed
display = triangles
border_color = black
color = none
line_width = 2
overlay_previous = share-y

[spacer]

[bed boundaries]
file = tad1_boundaries.bed
file_type = bed
height = 4
title = TAD Boundaries
min_value = 0
max_value = 2.5
```

Попробуем сделать по-другому (и немного другие параметры в ini, добавим трек CTCF):

```bash
hicFindTADs -m Kc167.corr2.20000.h5 --outPrefix tad2 --minDepth 60000 --maxDepth 200000  --step 20000 --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr
hicPlotTADs --tracks tracks2.ini --region chr2R:1000000-8000000 -o tads2.png
```

![tads2](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/tads2.png)

Здесь ТАДы размером от нескольких десятков килобаз, до 1.3 Мб. Также можно увидеть профиль ChIP-Seq для CTCF - многие пики приходятся на границы ТАД-ов. Теперь сравним ТАДы с этим профилем с помощью deeptools (e):

```bash
computeMatrix reference-point -S Kc167-CTCF.bigWig -R tad2_boundaries.bed --referencePoint center -a 200000 -b 200000 -out me.tab.gz -p 8
plotHeatmap -m me.tab.gz -out me.png --heatmapHeight 15 --heatmapWidth 10 --colorMap jet --sortRegions ascend --regionsLabel 'CTCF'
```

![me](/home/leo/ME/HSE/NGS/8/NGS_HW_Hi-C_data/me.png)

По оси Х здесь расстояние от границы ТАД-ов в мегабазах. По оси Y - регионы границ ТАД-ов, сортированные по среднему значению (снизу), сигнал CTCF (сверху). Паттерн распределения пиков CTCF вокруг границы ТАД-ов действительно такой, что наибольший сигнал приходится именно на саму границу, потом ослабляется при отходе от границы, а затем снова усиливается, но не достигает такого значения как на границы. 

То есть, основной вывод такой - на границе ТАД-ов плотность контактов уменьшается, что может свидетельствовать о обширных участках открытого хроматина, и это косвенно подтверждается данными ChIP-Seq для CTCF.