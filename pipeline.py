import argparse
import os
import subprocess
from helper_scripts.parseFastq import ParseFastQ

class ClinicalDataRow:

    def __init__(self, name, color, barcode):
        self.name = name
        self.color = color
        self.barcode = barcode

    def __repr__(self):
        return f"Name: {self.name} Barcode: {self.barcode} Color: {self.color}"


def main(fastq_file_name):
    # process_and_trim(fastq_file_name)
    # align()
    compress_alignments()


def compress_alignments():
    print("Compressing and sort alignments into BAM files")
    clinical_data = get_clinical_data()

    os.makedirs("bams", exist_ok=True)
    os.makedirs("sorted_bams", exist_ok=True)

    for clinical_data_row in clinical_data:
        with open(f"bams/{clinical_data_row.name}.bam", "w") as output_file:
            try:
                compress_command = ["samtools", "view", "-bS", f"sams/{clinical_data_row.name}.sam"]
                subprocess.run(compress_command, check=True, stdout=output_file)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
        sort_command = ["samtools", "sort", "-m", "100M", "-o", f"sorted_bams/{clinical_data_row.name}.sorted.bam", f"bams/{clinical_data_row.name}.bam"]
        subprocess.run(sort_command, check=True)
        index_command = ["samtools", "index", f"sorted_bams/{clinical_data_row.name}.sorted.bam"]
        subprocess.run(index_command, check=True)

def align():
    print("Aligning data...")
    clinical_data = get_clinical_data()
    initial_align_command = ["bwa", "index", "dgorgon_reference.fa"]
    try:
        subprocess.run(initial_align_command, check=True)
        print("Initial align command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        exit()

    os.makedirs("sams", exist_ok=True)

    for clinical_data_row in clinical_data:
        align_command = ["bwa", "mem", "dgorgon_reference.fa", f"fastqs/{clinical_data_row.name}_trimmed.fastq"]
        with open(f"sams/{clinical_data_row.name}.sam", "w") as output_file:
            try:
                subprocess.run(align_command, check=True, stdout=output_file)
                print("Align Command executed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                exit()


def process_and_trim(fastq_file_name):
    print("Processing and trimming data...")
    clinical_data = get_clinical_data()

    for clinical_data_row in clinical_data:
        print(f"Processing and trimming clinical data for {clinical_data_row.name}")
        fastq_file = ParseFastQ(fastq_file_name)

        directory = 'fastqs'
        file_path = os.path.join(directory, f"{clinical_data_row.name}_trimmed.fastq")

        os.makedirs(directory, exist_ok=True)
        fastq_write_file = open(file_path, 'w')

        for fastq_row in fastq_file:
            if(fastq_row.seq_barcode == clinical_data_row.barcode):
                fastq_write_file.write(fastq_row.seq_header + "\n")
                fastq_write_file.write(fastq_row.seq_trimmed + "\n")
                fastq_write_file.write(fastq_row.qual_header + "\n")
                fastq_write_file.write(fastq_row.qual_trimmed + "\n")


def get_clinical_data():
    clinical_data: list[ClinicalDataRow] = []
    try:
        clinical_data_file = open("harrington_clinical_data.txt", "r")
    except:
        print("ERROR: the clinical data file 'harrington_clinical_data.txt' must be in the same directory as this script")
        exit()

    for clinical_data_row in clinical_data_file.readlines()[1:]:
        clinical_data_row = clinical_data_row.split('\t')
        clinical_data.append(ClinicalDataRow(name=clinical_data_row[0].strip(), color=clinical_data_row[1].strip(), barcode=clinical_data_row[2].strip()))

    return clinical_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    args = parser.parse_args()
    main(args.fastq)
