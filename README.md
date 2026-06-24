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

      import csv
      import matplotlib.pyplot as plt
      from google.colab import files
      
      # ==============================================================================
      # STEP 1 & 2: Upload 3 File Berbeda & Masukkan ke dalam List
      # ==============================================================================
      list_sekuens = []
      jumlah_file = 3
      
      print(f"--- PROSES UPLOAD {jumlah_file} FILE FASTA ---")
      
      for i in range(1, jumlah_file + 1):
          print(f"\n[UPLOAD] Silahkan klik tombol di bawah untuk file ke-{i}:")
          uploaded = files.upload() # Memunculkan tombol upload
    
    for nama_file in uploaded.keys():
        print(f"-> Berhasil membaca: {nama_file}")
        
        # Decode data biner menjadi teks per baris
        isi_file = uploaded[nama_file].decode("utf-8").splitlines()
        
        header = ""
        sekuens_isi = []
        
        for baris in isi_file:
            baris = baris.strip()
            if not baris:
                continue
            
            if baris.startswith(">"):
                if header and sekuens_isi:
                    sekuens_penuh = "".join(sekuens_isi)
                    list_sekuens.append((header, sekuens_penuh))
                header = baris
                sekuens_isi = []
            else:
                sekuens_isi.append(baris.upper())
                
        if header and sekuens_isi:
            sekuens_penuh = "".join(sekuens_isi)
            list_sekuens.append((header, sekuens_penuh))

      print("\n" + "="*60)
      print(f"[SUKSES] Semua file terbaca. Total sekuens di List: {len(list_sekuens)}")
      print("="*60 + "\n")
      
      
      # ==============================================================================
      # STEP 3 & 4: Hitung Frekuensi (Dictionary) & Hitung GC Content
      # ==============================================================================
      hasil_analisis = []
      
      for header, seq in list_sekuens:
          # Menggunakan Dictionary untuk menghitung frekuensi (STEP 3)
          frekuensi = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
          for nukleotida in seq:
              if nukleotida in frekuensi:
                  frekuensi[nukleotida] += 1
            
    # Hitung GC Content
    total_g_c = frekuensi['G'] + frekuensi['C']
    total_basa = len(seq)
    gc_content = (total_g_c / total_basa) * 100 if total_basa > 0 else 0
    
    hasil_analisis.append({
        'Header': header,
        'Sekuens': seq,
        'Frekuensi': frekuensi,
        'GC_Content': gc_content
    })


      # ==============================================================================
      # STEP 4 & 5: Urutkan Berdasarkan GC Content & Tampilkan 3 Terbaik
      # ==============================================================================
      # Sorting descending berdasarkan nilai GC_Content
      hasil_terurut = sorted(hasil_analisis, key=lambda x: x['GC_Content'], reverse=True)
      
      print("=" * 60)
      print("           3 SEKUENS TERBAIK (GC CONTENT TERTINGGI)           ")
      print("=" * 60)
      
      top_3 = hasil_terurut[:3]
      for idx, data in enumerate(top_3, 1):
          print(f"Peringkat {idx}:")
          print(f"  Header    : {data['Header'][:80]}...") # Potong judul biar rapi di terminal
          print(f"  GC Content: {data['GC_Content']:.2f}%")
          print(f"  Frekuensi : {data['Frekuensi']}")
          print("-" * 60)
      
      
      # ==============================================================================
      # STEP 6: Visualisasi Grafik Batang (Bar Chart)
      # ==============================================================================
      # Mengambil ID pendek aksesi untuk label grafik (misal mengambil kata pertama setelah '>')
      nama_plot = [data['Header'].split()[0].replace('>', '') for data in hasil_terurut]
      nilai_gc = [data['GC_Content'] for data in hasil_terurut]
      
      plt.figure(figsize=(8, 5))
      bars = plt.bar(nama_plot, nilai_gc, color=['#008080', '#20b2aa', '#48d1cc'], edgecolor='black', width=0.5)
      
      plt.title('Perbandingan Nilai GC Content 3 Sekuens', fontsize=12, fontweight='bold')
      plt.xlabel('ID Aksesi Sekuens (NCBI)', fontsize=10)
      plt.ylabel('GC Content (%)', fontsize=10)
      plt.ylim(0, 100)
      
      for bar in bars:
          yval = bar.get_height()
          plt.text(bar.get_x() + bar.get_width()/2, yval + 2, f"{yval:.2f}%", ha='center', va='bottom', fontweight='bold')
      
      plt.tight_layout()
      plt.show()
      
      
      # ==============================================================================
      # STEP 7: Menuliskan Hasil Akhir ke File CSV
      # ==============================================================================
      nama_file_csv = "hasil_analisis_3_sekuens.csv"
      
      with open(nama_file_csv, mode='w', newline='', encoding='utf-8') as csv_file:
          fieldnames = ['Header', 'Panjang_Basa', 'A', 'T', 'C', 'G', 'GC_Content(%)']
          writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    
    writer.writeheader()
    for data in hasil_terurut:
        writer.writerow({
            'Header': data['Header'],
            'Panjang_Basa': len(data['Sekuens']),
            'A': data['Frekuensi']['A'],
            'T': data['Frekuensi']['T'],
            'C': data['Frekuensi']['C'],
            'G': data['Frekuensi']['G'],
            'GC_Content(%)': round(data['GC_Content'], 2)
        })

      print(f"\n[INFO] File output CSV berhasil dibuat: '{nama_file_csv}'")
