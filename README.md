# mini-project-bioinformatika-Rakan
# Simple Bioinformatics Pipeline: DNA Sequence Analysis & GC Content Calculation

## Deskripsi Proyek
Proyek ini adalah sebuah *pipeline* analisis sekuens bioinformatika sederhana yang dibangun menggunakan bahasa pemrograman **Python**. Program ini dirancang untuk melakukan pemrosesan massal (*batch processing*) terhadap file sekuens genetik berformat **FASTA** yang diunduh langsung dari database NCBI. 

Tujuan utama dari proyek ini adalah mengekstrak informasi nukleotida dari beberapa organisme/gen berbeda, menghitung frekuensi kemunculan masing-masing basa (A, T, C, G), menghitung persentase kandungan GC (*GC Content*), melakukan pengurutan (*sorting*) untuk menemukan sekuens dengan stabilitas termal tertinggi, serta menyajikan hasilnya dalam bentuk visualisasi grafik batang dan laporan data terstruktur format CSV.

---

## 🧬 Dasar Teori Singkat (Analisis Kontekstual)
Dalam studi bioinformatika dan genetika, **GC Content** (persentase pasangan basa Guanin dan Sitosin di dalam molekul DNA/RNA) merupakan parameter yang sangat krusial. 

Pasangan basa G-C dihubungkan oleh **3 ikatan hidrogen**, sedangkan pasangan A-T hanya dihubungkan oleh **2 ikatan hidrogen**. Akibatnya, sekuens DNA dengan kandungan GC yang tinggi memiliki stabilitas termal yang lebih kuat dan membutuhkan suhu hancur (*melting temperature* / $T_m$) yang lebih tinggi untuk memisahkan untai gandanya. Karakteristik ini sangat penting dalam optimasi desain primer PCR, studi adaptasi organisme termofilik, hingga anotasi genom.

---

## Alur Pipeline Program
Program berjalan secara berurutan (*pipeline*) melalui tahapan-tahapan berikut:



1. **Input & Parsing (Membaca File FASTA):** Program meminta *upload* file FASTA secara interaktif. Kode akan melakukan pembersihan data (*striping*) dan memisahkan string *header* (garis berawalan `>`) dari baris sekuens nukleotidanya.
2. **Data Storage (Penyimpanan List):** Seluruh pasangan data header dan sekuens utuh yang telah digabungkan disimpan ke dalam sebuah struktur data **List** di Python.
3. **Karakterisasi Nukleotida (Dictionary):** Menggunakan struktur data **Dictionary**, program melakukan iterasi tunggal pada string sekuens untuk menghitung frekuensi riil kemunculan basa Adenin (A), Timin (T), Sitosin (C), dan Guanin (G).
4. **Kalkulasi GC Content:** Menghitung rasio kandungan GC terhadap total panjang sekuens menggunakan rumus matematika sederhana:
   $$\text{GC Content (\%)} = \left( \frac{G + C}{A + T + C + G} \right) \times 100$$
5. **Sorting & Filter Top 3:** Seluruh data objek diurutkan secara *descending* (dari nilai terbesar ke terkecil) berdasarkan nilai GC Content untuk menyaring 3 sekuens terbaik teratas.
6. **Visualisasi (Matplotlib):** Membentuk grafik batang (*bar chart*) yang membandingkan kontras nilai persentase GC antar sekuens untuk memudahkan analisis visual.
7. **Export Data (CSV):** Menyimpan seluruh data tabel hasil kalkulasi matematis ke dalam file `hasil_analisis_3_sekuens.csv`.

---

## 📊 Dataset NCBI yang Digunakan
Proyek ini menguji 3 buah sekuens pendek yang merepresentasikan variasi genomik makhluk hidup:
1. **V00491.1** (*Homo sapiens* alpha 1 globin gene) - Representasi gen fungsional manusia.
2. **AH002844.2** (*Homo sapiens* insulin (INS) gene) - Representasi sekuens hormon manusia.
3. **NC_001422.1** (*Escherichia virus PhiX174*) - Representasi genom utuh dari virus bakteriofag pendek.

---

## Cara Menjalankan Program
1. Buka file notebook proyek (`.ipynb`) di Google Colab.
2. Jalankan cell kode program utama.
3. Saat program memunculkan tombol *Choose Files*, unggah file `.fasta` satu per satu sesuai instruksi di layar.
4. Program akan langsung mengeluarkan hasil *printout* peringkat di terminal, menampilkan grafik visual, dan membuat file CSV otomatis di direktori penyimpanan sementara Colab.

Note: maksimal 3 file FASTA
