# Analiz adımları: hücre hücre açıklama

Bu belge, Colab defterindeki her adımda ne yapıldığını, neden yapıldığını ve ne sonuç çıktığını anlatır. Amaç: **sağlıklı (kontrol) ve enfekte sığır meme dokusu arasında hangi genlerin farklı çalıştığını bulmak.**

## Temel fikir

RNA-seq, bir dokuda hangi genlerin ne kadar aktif olduğunu ölçer. Aktif bir gen çok sayıda mRNA kopyası üretir. Bu mRNA'lar parçalanıp okunur ve bir gene ait okuma sayısı o genin aktivitesinin göstergesi olur. İki grupta bu sayıları istatistiksel olarak karşılaştırırız.

```
FASTQ (ham okumalar) -> kalite kontrol -> temizleme -> gen sayıları -> DESeq2 -> sonuçlar
```

## Küçük sözlük

| Terim | Anlamı |
|---|---|
| FASTQ | Ham okuma dosyası. Her okuma 4 satırdır: ad, baz dizisi (A, C, G, T), `+`, her baz için kalite karakteri |
| Paired-end | Her RNA parçasının iki ucu okunur, bu yüzden her örnek için iki dosya var (`_1`, `_2`) |
| Kalite skoru (Phred) | Bazın doğru okunma güvenilirliği. `J` çok yüksek, `#` çok düşük (Q2) |
| Adaptör | Laboratuvarda parçanın ucuna eklenen sabit dizi. Parça kısaysa okuma adaptöre taşar, bu kısım gen bilgisi taşımaz |
| Sayım tablosu (count matrix) | Satırlar genler, sütunlar örnekler, hücreler o gene düşen okuma sayısı |
| log2 kat değişimi (log2FC) | 1 = 2 kat artış, -1 = yarıya düşüş, 2 = 4 kat artış |
| padj | 14.000 gen aynı anda test edildiği için düzeltilmiş p değeri. 0,05'in altı anlamlı kabul edilir |

---

## Bölüm 1: Hazırlık

### 1.1 Drive'ı bağla, klasörleri hazırla
Colab'ın geçici diski oturum kapanınca silinir. Sonuçları kaybetmemek için Google Drive'ı bağlıyoruz. `BASE` kalıcı klasörümüz (Drive), `RAW` ham veri için geçici klasör.
**Sonuç:** `hazır`

### 1.2 Örnek listesini oku
`samples.csv` dosyasını okur. İlk plandaki 8 örnek (4 inek x enfekte/kontrol) listelenir, SRA dosya boyutları 3,2 ile 4,8 GB arasında.
**Not:** DESeq2 aşamasında 6 ineğin 12 örneği kullanıldı (sayım tablosu hazır olduğu için ek maliyet yoktu).

### 1.3 Araçları kur
`apt-get` ile dört program kurulur: **sra-toolkit** (veri indirme), **pigz** (hızlı sıkıştırma), **FastQC** (kalite kontrol), **Trimmomatic** (temizleme). Colab her yeni oturumda sıfırlandığı için bu hücre her seferinde çalıştırılır.
**Sonuç:** `fasterq-dump 3.0.3`, diskte 206 GB boş alan.

---

## Bölüm 2: Veri indirme (test örneği SRR2960032)

### 1.4 Tek örnek testi
- `prefetch`: Örneği NCBI'nin SRA arşivinden indirir (sıkıştırılmış `.sra` formatı).
- `fasterq-dump --split-files`: `.sra` dosyasını iki FASTQ dosyasına çevirir.
- `pigz`: FASTQ'ları sıkıştırır (diskte yer kazanmak için).

**Sonuç:** 25.250.142 okuma çifti (50.500.284 okuma), iki dosya, her biri 2,0 GB.

### 1.5 Kontrol
Dosyanın ilk okumasını gösterir ve satır sayısını sayar.
**Sonuç:** İlk okuma 101 baz uzunluğunda. Satır sayısı 101.000.568, bu 4 x 25.250.142 eder, yani dosya eksiksiz ve bozulmamış.

---

## Bölüm 3: Kalite kontrol

### FastQC (ham veri)
FastQC, FASTQ dosyalarını 10 farklı açıdan inceler ve her biri için PASS (geçti), WARN (uyarı) ya da FAIL (kaldı) verir. Her okuma dosyası için bir HTML raporu üretir.

### Zip içinden adaptör ve aşırı temsil edilen dizileri okuma
FastQC'nin metin çıktısından iki bölümü yazdırır.
**Sonuç (önemli bulgu):**
- Okuma 1'de en sık geçen üç dizi **TruSeq adaptör** ve PCR primeri: toplam okumaların yaklaşık **%9'u**. Bu okumalar gen bilgisi taşımıyor.
- Okumaların sonuna doğru adaptör oranı: okuma 1'de yaklaşık %5, okuma 2'de yaklaşık %3.
- Kalite, GC oranı (%49) ve okuma uzunluğu sağlıklı. Düşük kaliteli dizi yok.
- Sequence Duplication FAIL: RNA-seq'te çok ifade edilen genler çok tekrar üretir, bu normaldir.

### Trimmomatic (temizleme)
Java ile çalıştırılır. Ayarlar:

| Ayar | Ne yapar |
|---|---|
| `ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True` | Adaptör dizilerini bulup keser |
| `LEADING:3`, `TRAILING:3` | Okumanın başındaki ve sonundaki çok düşük kaliteli bazları atar |
| `SLIDINGWINDOW:4:15` | 4 bazlık pencerede ortalama kalite 15'in altına inerse okumanın kalanını keser |
| `MINLEN:36` | Kırpıldıktan sonra 36 bazdan kısa kalan okumaları atar |

**Sonuç (17,2 dakika):** 25.250.142 çiftin %82,73'ü iki ucuyla sağlam kaldı, %4,43'ü yalnızca ileri, %1,08'i yalnızca geri okuma olarak kaldı, %11,76'sı elendi. Elenenlerin büyük kısmı zaten sadece adaptörden oluşan okumalardı.

### FastQC (temizlenmiş veri)
Aynı kontrol temizlenmiş dosyalarda tekrarlanır.
**Sonuç:** **Adapter Content iki okumada da PASS.** Temizlik işe yaradı. Per base sequence content FAIL ve GC WARN gibi bazı uyarılar RNA-seq kütüphanelerinde sık görülür ve gen düzeyi analizi genelde bozmaz.

---

## Bölüm 4: Diferansiyel ifade analizi

Not: Bu bölüm, test örneğinden değil, GEO'da yazarların yayımladığı **hazır sayım tablosundan** (`GSE75379_count.txt.gz`) yapıldı. Test örneği yalnızca indirme ve kalite kontrol adımlarını göstermek içindi.

### Dosyaları bul ve girdileri hazırla (`01_prepare_inputs.py`)
- Sayım tablosunu (24.616 gen x 36 örnek, meme ve karaciğer) okur.
- Metadata dosyasına göre yalnızca bizim 12 meme örneğini (24. saat) seçer.
- **Sayıların tam sayı olup olmadığını kontrol eder.** DESeq2 ham sayım ister, normalize değerler (örneğin RPKM) uygun değildir. Kontrol geçti, yani tablo ham sayım.
- `counts.csv` ve `coldata.csv` dosyalarını yazar.

**Sonuç:** 12 örnek, 6 inek. Örnek başına toplam okuma 11,5 ile 32,1 milyon arasında. Bu fark normal ve DESeq2 bunu dengeler.

### DESeq2 (`02_deseq2.R`)
Adım adım:

1. **Filtre:** En az 3 örnekte 10 okumadan az olan genler atılır. Kalan: **14.369 gen**.
2. **Model:** `~ inek + durum`. Aynı ineğin enfekte ve kontrol çeyreği eşleştirilir, böylece ineğe özgü farklar karşılaştırmayı bozmaz.
3. **`estimating size factors`:** Örnekler arasındaki toplam okuma farkını dengeler (normalizasyon).
4. **`estimating dispersions`:** Her genin örnekler arasındaki doğal değişkenliğini tahmin eder. Fark bundan büyükse anlamlıdır.
5. **`fitting model and testing`:** Her gen için enfekte ve kontrol arasında test yapar, p değerini hesaplar ve padj'ye düzeltir.

**Sonuç:**
- padj < 0,05: **3.826 gen** (2.079 artan, 1.747 azalan).
- Buna ek olarak |log2FC| ≥ 1 şartıyla: **1.931 gen** (1.236 artan, 695 azalan).

**Grafikler:**
- **PCA:** Örnekleri tüm genlere bakarak 2 boyuta indirger. Birinci bileşen (%67) kontrol ve enfekte örnekleri tamamen ayırıyor. Bu, farkın gerçek olduğunun en güçlü işareti.
- **Volcano:** Her nokta bir gen. Sağ taraf artanlar, sol taraf azalanlar, yukarısı daha anlamlı.
- **Isı haritası:** En anlamlı 30 genin örneklerdeki ifadesi. Örnekler iki gruba ayrılıyor ve inekler gruplara karışmış.

### Gen adlarını ekleme ve en önemli genler
Ensembl kimliklerini (ENSBTAG...) Ensembl REST ile gen adlarına çevirir ve tabloyu kaydeder.

**Sonuç (özet):**

| | Genler | Süreç |
|---|---|---|
| Artan | CXCR1, CSF3R, CCR1, S100A9, SAA2, HP, CD68, CD163, CFB, IL1RAP, SLC11A1 | Nötrofil ve makrofaj göçü, akut faz yanıtı, kompleman sistemi |
| Azalan | FASN, LPL, BTN1A1, BTN2A2, THRSP | Süt yağı sentezi ve lipit metabolizması |

Bazı genlerin adı bulunamadı (eski kimlikler), onlar kimlikle bırakıldı.

---

## Bölüm 5: Çıkmaz sokaklar (defterden silinebilir)

Defterde sonuca katkısı olmayan hücreler de var, bunlar süreçte yaşanan sorunlardan kalma:

- **`NameError: name 'QC' is not defined`:** Colab oturumu yenilenince değişkenler silinir.
- **`trimmomatic: command not found`:** Trimmomatic kurulu ama bu adla bir komut olarak tanımlı değildi. Çözüm: `java -jar` ile doğrudan çalıştırmak.
- **Salmon hücreleri (kurulum, referans indirme, indeks, `quant`, hata ayıklama):** Ham okumalardan kendi gen sayımımızı çıkarmaya çalıştık. Colab oturumları sıfırlandığı için tamamlanmadı. Bunun yerine yazarların sayım tablosunu kullandık.
- Boş hücreler.

## Bu projeden çıkan dersler

- Colab'ın disk ve bellek durumu oturum kapanınca silinir. Önemli çıktıları hep Drive'a kaydet.
- Ham veriyle çalışmak saatler sürer. İlk deneme için bir örnek üzerinde hattı test etmek akıllıca.
- Kalite kontrol bir kez çalıştırıp geçilen bir adım değil: temizledikten sonra tekrar bakarak işe yaradığını doğruladık.
- PCA gibi basit bir grafik, 14.000 genlik bir sonucun güvenilirliğini bir bakışta gösterebilir.
