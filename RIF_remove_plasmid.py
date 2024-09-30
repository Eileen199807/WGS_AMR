# Usage: python plasmid.py
import os
import subprocess
import sys
import shutil
import pandas as pd
import glob
from Bio import SeqIO
import re
import shlex

# .bed files generation, used for replace the plasmids sequences with NNNNNNN
def generate_bed_files(plasmid_output_folder):
    for folder in os.listdir(plasmid_output_folder):
        folder_path = os.path.join(plasmid_output_folder, folder)
        if os.path.isdir(folder_path):
            tsv_files = glob.glob(os.path.join(folder_path, "results_tab.tsv"))
            for tsv_file in tsv_files:
                data = pd.read_csv(tsv_file, sep="\t")
                Position = data['Position in contig'].tolist()
                Contig = data['Contig'].tolist()
                result_list = []
                for contig, position in zip(Contig, Position):
                    contig_id = contig if isinstance(contig, int) else contig.split()[0]
                    start, end = position.split('..')
                    start = int(start) - 1
                    result = f"{contig_id} {start} {end}"
                    result_list.append(result)
                data = [result.split() for result in result_list]
                df = pd.DataFrame(data)
                output_file = os.path.join(folder_path, f"{folder}.bed")
                df.to_csv(output_file, sep='\t', header=False, index=False)
                print(f"Saved {output_file}")

# Count bases in each contig for FASTA files and calculate of ther percentage of plasmids in each contig
def process_results(plasmid_output_folder, fasta_folder):
    def get_contig_base_count(fasta_file, target_contig_id):
        base_count = 0
        with open(fasta_file, "r") as file:
            for record in SeqIO.parse(file, "fasta"):
                contig_id = record.id
                if contig_id == target_contig_id:
                    base_count = len(record.seq)
                    break
        return base_count
    def extract_first_part(contig_info):
        if isinstance(contig_info, int):
            return contig_info
        else:
            return contig_info.split(" ")[0]
    os.makedirs("Plasmidfinder_process", exist_ok=True)
    for folder in os.listdir(plasmid_output_folder):
        folder_path = os.path.join(plasmid_output_folder, folder)
        if os.path.isdir(folder_path):
            tsv_file = os.path.join(folder_path, "results_tab.tsv")
            fasta_file = os.path.join(fasta_folder, f"{folder}.fasta")
            data = pd.read_csv(tsv_file, sep="\t")
            data['Contig ID'] = data['Contig'].apply(extract_first_part)
            data['Base count'] = data['Contig ID'].apply(lambda x: get_contig_base_count(fasta_file, str(x)))
            data['Plasmid length'] = data['Position in contig'].apply(lambda x: int(x.split('..')[1]) - int(x.split('..')[0]) + 1)
            data['Plasmid length/Base count'] = data['Plasmid length'] / data['Base count']           
            output_file = os.path.join("Plasmidfinder_process", f"{folder}_results.csv")
            data.to_csv(output_file, index=False)
            print(f"Saved {output_file}")

def combine_and_calculate(processed_folder):
    files = os.listdir(processed_folder)
    csv_files = [f for f in files if f.endswith('.csv')]
    data_frames = [] 
    for csv_file in csv_files:
        data = pd.read_csv(os.path.join(processed_folder, csv_file))
        data['Filename'] = csv_file.rsplit('_', 1)[0]
        data_frames.append(data)
        combined_data = pd.concat(data_frames, ignore_index=True)
        result = combined_data.groupby(['Contig ID', 'Filename'])['Plasmid length/Base count'].sum().reset_index()
        result_data = result[['Contig ID', 'Plasmid length/Base count', 'Filename']]
        Delete_contig = result_data[result_data['Plasmid length/Base count'] > 0.2] 
    Delete_contig.to_csv("Contig_deletion.csv", index=False)    
    result_data.to_csv("Combined_plasmid_result.csv", index=False)

def is_contig_to_remove(seq_id, contig_number):
    pattern = f"^{contig_number}$"
    return re.search(pattern, seq_id) is not None

def remove_contigs(input_path, table_data):
    for index, row in table_data.iterrows():
        input_file = os.path.join(input_path, row["Filename"] + ".fasta")
        contig_to_remove_number = row["Contig ID"]
        print(contig_to_remove_number)
        if os.path.exists(input_file):
            sequences = list(SeqIO.parse(input_file, "fasta"))
            filtered_sequences = [seq for seq in sequences if not is_contig_to_remove(seq.id, contig_to_remove_number)]
            SeqIO.write(filtered_sequences, input_file, "fasta")
        else:
            print(f"File {input_file} not found.")

def mask_fasta_with_bedtools(fasta_file, bed_file, output_file):
    bedtools_cmd = f"/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/Softwares/bedtools2/bin/bedtools maskfasta -fi {fasta_file} -bed {bed_file} -fo {output_file}"
    subprocess.run(bedtools_cmd, shell=True, check=True)

def process_all_fna_files(process_input_folder, bed_folder):
    for file_name in os.listdir(process_input_folder):
        if file_name.endswith(".fasta"):
            sample_id = file_name[:-6]
            fasta_file = os.path.join(process_input_folder, file_name)
            bed_file = os.path.join(plasmid_output_folder, sample_id, f"{sample_id}.bed")
            output_file = os.path.join(process_input_folder, f"{sample_id}_noplasmid.fasta")
            if os.path.exists(bed_file):
                mask_fasta_with_bedtools(fasta_file, bed_file, output_file)
                os.remove(fasta_file) 
            else:
                print(f"Bed file {bed_file} not found.")
                os.rename(fasta_file, output_file)  

if __name__ == "__main__":
    
    plasmid_output_folder="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/TEST/SAURS_Plasmidfinder"
    folder_path="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/TEST/SAURS_FASTA"
    
    generate_bed_files(plasmid_output_folder)
    process_results(plasmid_output_folder, folder_path)
    
    processed_folder = "/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/TEST/Plasmidfinder_process"
    combine_and_calculate(processed_folder)
    
    subprocess.run(["cp", "-r", "SAURS_FASTA", "PROCESS_FASTA"])
    process_path = os.path.abspath("PROCESS_FASTA")

    table_data = pd.read_csv("Contig_deletion.csv")
    remove_contigs(process_path, table_data)
    process_all_fna_files(process_path, plasmid_output_folder)

   