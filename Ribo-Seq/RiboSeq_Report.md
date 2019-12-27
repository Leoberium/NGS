# RiboSeq

### Лев Мазаев | мАДБМ18

Считал все на локальной машине.

### Этап 1: скачивание и распаковка файлов образцов

Первая задача - скачать SRA-файлы. Я решил выбрать образцы с нокдауном метилтрансферазы METTL3. На  GEO-странице эксперимента внизу можно перейти в SRA Run Selector:

![src1](/home/leo/ME/HSE/NGS/5/src1.png)

Там мы видим все образцы (Samples) и их соответствующие SRR-записи:

![src2](/home/leo/ME/HSE/NGS/5/src2.png)

Скачиваем SRR586579(4|5) (GSM271736(0|1) - Ribo_seq_shMETTL3_rep(1|2)) и SRR586580(1|2) (GSM271736(7|8) - RNA_seq_shMETTL3_rep(1|2)):

```bash
mkdir ~/BioData/RiboSeq
cd ~/BioData/Riboseq
prefetch -v SRR5865794
prefetch -v SRR5865795
prefetch -v SRR5865801
prefetch -v SRR5865802
```

Далее, в процессе скачивания посмотрим, что делают команды `esearch` и `efetch`:

![src3](/home/leo/ME/HSE/NGS/5/src3.png)

Из них при помощи ruby-скрипта можно получить табличку для выяснения номеров SRA, но не стал этого делать, так как всё есть в SRA Run Selector. 

Распакуем скачанные SRA-архивы (предварительно собрав их в одной папке) с помощью `parallel-fastq-dump` (работает гораздо быстрее, чем обычный `fastq-dump`):

```bash
parallel-fastq-dump -s *.sra -t 8 --split-files
```

Сравнение размеров файлов:

![src4](/home/leo/ME/HSE/NGS/5/src4.png)

### Этап 2: отрезка адаптеров и QC

Первичный анализ последовательностей:

```bash
fastqc -t 8 ~/BioData/RiboSeq/*.fastq
```

Длина рида в первых репликах - 49, во вторых - 51. Заметен полиА-адаптер в первых репликах (рис. для RiboSeq rep1):

![fc1](/home/leo/ME/HSE/NGS/5/fc1.png)

Распределение длин прочтений (RiboSeq rep1):

![fc2](/home/leo/ME/HSE/NGS/5/fc2.png)

RNASeq rep1:

![fc3](/home/leo/ME/HSE/NGS/5/fc3.png)

Для вторых реплик тоже самое с пиком на 51. Теперь произведём отрезку адаптеров.

![src5](/home/leo/ME/HSE/NGS/5/src5.png)

Адаптер первой реплики - АААААААААА:

```bash
cutadapt -a AAAAAAAAAA -q 20 --minimum-length 20 -j 8 -o riboseq1.trimmed.fastq.gz --trimmed-only SRR5865794_1.fastq
cutadapt -a AAAAAAAAAA -q 20 --minimum-length 20 -j 8 -o rnaseq1.trimmed.fastq.gz SRR5865801_1.fastq
```

Адаптер второй реплики - CTGTAGGCACCATCAAT (minimum length 30, так как потом будем ещё отрезать баркоды):

```bash
cutadapt -a CTGTAGGCACCATCAAT -q 20 --minimum-length 30 -j 8 -o riboseq2.trimmed.fastq.gz --trimmed-only SRR5865795_1.fastq
cutadapt -a CTGTAGGCACCATCAAT -q 20 --minimum-length 30 -j 8 -o rnaseq2.trimmed.fastq.gz SRR5865802_1.fastq
```

Далее, нужно дедуплицировать вторые реплики:

```bash
seqkit rmdup -s -j 8 riboseq2.trimmed.fastq.gz | gzip -c > riboseq2.dedup.fastq.gz
```

[INFO] 18088740 duplicated records removed. То есть 46% записей было удалено из исходных 38801308 в `riboseq2.trimmed.fastq.gz`.

```bash
seqkit rmdup -s -j 8 rnaseq2.trimmed.fastq.gz | gzip -c > rnaseq2.dedup.fastq.gz
```

[INFO] 9371952 duplicated records removed. То есть 36% записей было удалено из исходных 25769506 в `rnaseq2.trimmed.fastq.gz`

Теперь удалим баркоды:

```bash
cutadapt -u 6 -j 8 riboseq2.dedup.fastq.gz | cutadapt -u -4 -j 8 - -o riboseq2.final.fastq.gz
cutadapt -u 6 -j 8 rnaseq2.dedup.fastq.gz | cutadapt -u -4 -j 8 - -o rnaseq2.final.fastq.gz
```

FastQC после всех операций. RiboSeq rep1:

![fc4](/home/leo/ME/HSE/NGS/5/fc4.png)

RNASeq rep1:

![fc5](/home/leo/ME/HSE/NGS/5/fc5.png)

RiboSeq rep2:

![fc6](/home/leo/ME/HSE/NGS/5/fc6.png)

RNASeq rep2:

![fc7](/home/leo/ME/HSE/NGS/5/fc7.png)

Видим, что характерная длина рибосомных футпринтов - 25-32 пн. Из первой реплики RNASeq можно сделать вывод о том, что авторы старались сделать фрагменты РНК близкими по длине к рибосомным.

### Этап 3: проверим долю рибосомной РНК

Построим индекс `bowtie` (скачал rRNA_euk.fasta с сервера):

```bash
bowtie-build --threads 8 rRNA_euk.fasta rRNA
```

Получим оценки доли рРНК для каждого из образцов, полученных выше:

![src6](/home/leo/ME/HSE/NGS/5/src6.png)

### Этап 4: картирование на реальный геном

Использовал готовые файлы.

### Этап 5: фазирование прочтений и метагенные профили

Фазирование по RiboSeq нокдауну METTL3 без aggregate:

```bash
psite ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt psite_test --countfile_format BAM --count_files ~/BioData/RiboSeq/olbams/METTL3_ribo_Coots2017_m_r1.bam ~/BioData/RiboSeq/olbams/METTL3_ribo_Coots2017_m_r2.bam --min_length 20 --max_length 40 --constrain 10 18 --min_count 10 --default 14
```



![wo_aggregate](/home/leo/ME/HSE/NGS/5/psite/wo_aggregate.png)

Фазирование по RNASeq (от 30 до 50):

```bash
psite ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt psite_test --countfile_format BAM --count_files ~/BioData/RiboSeq/olbams/METTL3_rna_Coots2017_m_r1.bam ~/BioData/RiboSeq/olbams/METTL3_rna_Coots2017_m_r2.bam --min_length 30 --max_length 50 --aggregate --constrain 10 18 --min_count 10 --default 14
```



![rnaseq](/home/leo/ME/HSE/NGS/5/psite/rnaseq.png)

Фазирование по всем RiboSeq-образцам:

```bash
psite ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt psite_test --countfile_format BAM --countfiles $(find ~/BioData/RiboSeq/bams/ -name '*ribo*.bam' | tr '\n' ' ') --min_length 20 --max_length 40 --aggregate --constrain 10 18 --min_count 10 --default 14
```



![psite_test_p_offsets](/home/leo/ME/HSE/NGS/5/psite/psite_test_p_offsets.png)

Проверка фазирования:

```bash
phase_by_size --countfile_format BAM --count_files $(find ~/BioData/RiboSeq/bams/ -name '*ribo*.bam' | tr '\n' ' ') --fiveprime_variable --offset ./psite/psite_test_p_offsets.txt --min_length 25 --max_length 31 ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt phas_by_size_test
```

![phas_by_size_test_phasing](/home/leo/ME/HSE/NGS/5/phas_by_size_test_phasing.png)

Получилось достаточно правильное, большая часть ридов ложится на frame0.

Нефазированный метагенный профиль METTL3 rep2 (длина 25-31):

```bash
metagene count --countfile_format BAM --count_files ~/BioData/RiboSeq/bams/METTL3_ribo_Coots2017_m_r2.bam --fiveprime --min_length 25 --max_length 31 --min_count 10 --use_mean --landmark Start ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt metagene_counts_test
```

![METTL3_rep2_unphased](/home/leo/ME/HSE/NGS/5/METTL3_rep2_unphased.png)

Фазированный метагенный профиль того же образца:

```bash
metagene count --countfile_format BAM --count_files ~/BioData/RiboSeq/bams/METTL3_ribo_Coots2017_m_r2.bam --fiveprime_variable --offset ./psite/psite_test_p_offsets.txt --min_length 25 --max_length 31 --min_count 5 --use_mean --landmark Start ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt metagene_counts_test
```

![METTL3_rep2_phased](/home/leo/ME/HSE/NGS/5/METTL3_rep2_phased.png)Метагенный профиль RNASeq METTL3 rep2:

```bash
metagene count --countfile_format BAM --count_files ~/BioData/RiboSeq/bams/METTL3_rna_Coots2017_m_r2.bam --fiveprime --min_length 35 --max_length 45 --min_count 10 --use_mean --landmark Start ~/BioData/RiboSeq/plastidmetagene/mouse_start_rois.txt metagene_counts_test
```



![METTL3_rna](/home/leo/ME/HSE/NGS/5/METTL3_rna.png)

### Этап 6: получение bedGraph и визуализация

Сделаем профили RiboSeq и RNASeq для METTL3 rep2:

```bash
make_wiggle -o METTL3_ribo2 --count_files ~/BioData/RiboSeq/bams/METTL3_ribo_Coots2017_m_r2.bam --normalize --min_length 25 --max_length 31 --fiveprime_variable --offset ../psite/psite_test_p_offsets.txt 
make_wiggle -o METTL3_rna2 --count_files ~/BioData/RiboSeq/bams/METTL3_rna_Coots2017_m_r2.bam --normalize --center
bedtools unionbedg -i METTL3_ribo2_fw.wig METTL3_ribo2_rc.wig | csvtk mutate2 -H -t -L 5 -e '$4+$5' | cut -f 4-5 --complement > ribo.bedGraph
bedtools unionbedg -i METTL3_rna2_fw.wig METTL3_rna2_rc.wig | csvtk mutate2 -H -t -L 5 -e '$4+$5' | cut -f 4-5 --complement > rna.bedGraph
```

Визуализируем ген *ATF4* с помощью `svist4get`:

```bash
svist4get -gtf ~/BioData/Annotations/gencode.vM23.basic.annotation.gtf -fa ~/BioData/GRCm38.primary_assembly.genome.fa -bg ribo.bedGraph rna.bedGraph -g ENSMUSG00000042406.8
```



![svist4get](/home/leo/ME/HSE/NGS/5/svist4get.png)

Как можно заметить, данный ген имеет преждевременную рамку считывания в 5'-UTR.